#Takes a directory of cmsearch outputs, one for each calibrated model.  
#
#Outputs a dataframe of scores where the rows are query sequences and column names are
#sequences that have a match with that query.

#Usage: python parsecmsearch.py <directory of cmsearchoutputs> <outfile.txt>


import sys
import os
import numpy as np

def parseCMsearch(cmsearchoutput):
    #Parse output from cmsearch and return dictionary of hits and scores
    hitdict = {} #{sequencehit : [rank, evalue, score]}
    scorelines = []
    querydict = {} #{query : {sequencehit : [rank, evalue, score], sequencehit2 : [rank, evalue, score]}}
    infh = open(cmsearchoutput, 'r')
    for line in infh:
        #Only deal with non-blank lines
        if line not in ['\n', '\r\n']:
            line = line.strip().split()
            #Get query line
            if len(line) > 0 and line[0] == 'Query:':
                query = line[1].split(';')[0]
                querydict[query] = []
            #Get rank/sequence lines
            if len(line) == 13:
                if line[0][0] == '(':
                    if line[10] == '+' or '3\'' or '5\'':
                        if line[1] == '!' or '?':
                            scorelines.append(line)

    for line in scorelines:
        rank = int(line[0].replace('(','').replace(')',''))
        evalue = float(line[2])
        score = float(line[3])
        sequencehit = line[5]
        #If two hits in the same sequence, only take the higher scoring one
        if sequencehit not in hitdict:
            hitdict[sequencehit] = [rank, evalue, score]

    querydict[query] = hitdict
    return querydict

def compilequerydicts(cmsearchoutputdirectory):
    querydicts = []
    #Run parseCMsearch on all files in cmsearchoutputdirectory
    cmsearchoutputdirectory = os.path.abspath(cmsearchoutputdirectory)
    cmsearches = [os.path.join(cmsearchoutputdirectory, cmsearch) for cmsearch in os.listdir(cmsearchoutputdirectory) if os.path.isfile(os.path.join(cmsearchoutputdirectory, cmsearch))]

    for cmsearch in cmsearches:
        querydict = parseCMsearch(cmsearch)
        querydicts.append(querydict)

    #[{query1 : {sequencehit : [rank, evalue, score], sequencehit2 : [rank, evalue, score]}}, {query2 : {seqhit1 : [rank, evalue,
    #score], sequencehit2 : [rank, evalue, score]}}]
    return querydicts

def makecmscorematrix(querydicts, outfile):
    #Takes list of querydicts as input.  Outputs dataframe of scores where rows are queries and columns are searchhits
    
    #Querydicts is list of {query: {searchhit : [rank, evalue, score]}} dictionaries
    #Queries are in query.keys()[0].  They get added to [queries] in order that they appear in querydicts.
    #Get all searchhits...these are keys of the innermost dictionaries...equivalent to seqhit1, seqhit2 above
    queries = []
    searchhits = []
    for query in querydicts:
        queries.append(query.keys()[0])
        for value in query.values():
            for searchhit in value.keys():
                if searchhit not in searchhits:
                    searchhits.append(searchhit)

    #Turn querydicts into an array.  Colnames will be search hits in the order they are in [searchhits].
    #Row names will be queries in the order they are in [querydicts].
    #Each query must also be in searchhits.  This will happen if you a model's own sequence is
    #included in the sequences to search through. It should have a VERY high score.
    list_of_lists = []
    for querydict in querydicts:
        for query_name, query in sorted(querydict.items()):
            scores = []
            for searchhit in searchhits:
                try:
                    scores.append(query[searchhit][2])
                except KeyError:
                    scores.append(0)
                    
        list_of_lists.append(scores)

    querydicts_array = np.array(list_of_lists)

    #Output array to file
    tempoutfile = 'temp_' + outfile
    np.savetxt(tempoutfile, querydicts_array, delimiter = '\t', fmt = '%2f')

    #Put rownames(queries) in order as first column
    rows = [line.split('\t') for line in file(tempoutfile)]
    cols = zip(*rows)
    cols.insert(0, queries)
    rows = zip(*cols)
    file(outfile, 'w').writelines(['\t'.join(row) for row in rows])

    #Put colnames(searchhits) in order at top of output
    #Need blank item in searchhits because first column is now rownames!
    searchhits.insert(0,'')
    with open(outfile, 'r+') as f:
        old = f.read()
        f.seek(0)
        f.write(('\t').join(searchhits) + '\n' + old)

    os.remove(tempoutfile)
    
        

if __name__ == '__main__':
    querydicts = compilequerydicts(sys.argv[1])
    makecmscorematrix(querydicts, sys.argv[2])
