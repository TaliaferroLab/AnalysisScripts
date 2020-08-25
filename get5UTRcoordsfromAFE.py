#Takes a MISO gff of AFE events and a reference gff.  Reference must have start codons explicitly annotated as such.  Searches for start codons in exons of each AFE.  Starts with most gene-proximal exons.  If none, moves to the next outer exon.  If >1 in an exon, takes the most gene-proximal start to give the longest UTR.  If none in the entire AFE, prints a warning and skips.  Returns gff of 5' UTR portions of each ALE.  Child exons of each UTR are truncated so that they do not exceed boundaries of 5' UTR.

#Usage python get5UTRcoords.py <AFE.gff> <reference.gff> <outfile.gff>

import sys
import os
from operator import itemgetter
import gffutils

gff_fn = sys.argv[1]
db_fn = os.path.basename(gff_fn) + '.db'

if os.path.isfile(db_fn) == False : #if database doesn't exist, create it
    gffutils.create_db(gff_fn, db_fn)

db = gffutils.FeatureDB(db_fn)

def get5UTRcoords(AFEgff, startcodongff):

    startcodonfh = open(startcodongff, 'r')
    AFEfh = open(AFEgff, 'r')
    
    startcodons = []
    currentisoform = []
    AFEisoforms = {} #{(chrm, strand) : [exon1start, exon1stop], [exonNstart, exonNstop]}
    AFEs = []
    startlist = []
    AFEswithStarts = 0
    UTRgffAll = []
   
    for line in startcodonfh:
        line = line.strip()
        line = line.split('\t')
        if line[2] == 'start_codon':
            startcodons.append(line)

    for line in AFEfh:
        line = line.strip()
        line = line.split('\t')
        AFEs.append(line)

    startcodonfh.close()
    AFEfh.close()
            
    #Make dictionary {(chrm,strand) : [exon1start, exon1stop], ... , [exonNstart, exonNstop]}
    for idx, line in enumerate(AFEs):
        chrm = line[0]
        cat = line[2]
        start = int(line[3])
        stop = int(line[4])
        strand = line[6]
        if cat == 'mRNA':
            ID = line[8]

        if cat == 'mRNA':
            isoformstart = start
            isoformstop = stop

        if cat == 'exon':
            
            currentisoform.append([start, stop])
            if idx + 1 == len(AFEs): #if at the end of the list
                if strand == '+':
                    currentisoform = sorted(currentisoform, key = itemgetter(0), reverse = True) #exon1 is closest to CDS to give longest UTR
                elif strand == '-':
                    currentisoform = sorted(currentisoform, key = itemgetter(0))
                AFEisoforms[tuple([chrm, isoformstart, isoformstop, strand, ID])] = currentisoform
                
               
            elif idx + 1 < len(AFEs):
                nextcat = AFEs[idx + 1][2]
                if nextcat != 'exon': #if next entry is not an exon entry
                    if strand == '+':
                        currentisoform = sorted(currentisoform, key = itemgetter(0), reverse = True)
                    elif strand =='-':
                        currentisoform = sorted(currentisoform, key = itemgetter(0))
                    AFEisoforms[tuple([chrm, isoformstart, isoformstop, strand, ID])] = currentisoform
                    currentisoform = []
    
    
    for AFEisoform in AFEisoforms:
        d = {}
        startlist = []
        AFEchrm = AFEisoform[0]
        AFEstart = AFEisoform[1]
        AFEstop = AFEisoform[2]
        AFEstrand = AFEisoform[3]
        AFEID = AFEisoform[4]
        numExons = len(AFEisoforms[AFEisoform])
        for idx, exon in enumerate(AFEisoforms[AFEisoform]): #dictionary {'exonstartN' : 12345, 'exonstopN': 23456}
            d['exonstart%s' % str(idx+1)] = exon[0]
            d['exonstop%s' % str(idx+1)] = exon[1]

       
        

        for startcodon in startcodons:
            startcodonchrm = startcodon[0]
            startcodonstart = int(startcodon[3])
            startcodonstop = int(startcodon[4])
            startcodonstrand = startcodon[6]
            currentexon = 1 #start with exon 1.  This is the one most proximal to the gene to give the longest UTRs.
            

            if AFEchrm == startcodonchrm and AFEstrand == startcodonstrand:
                for exon in d:
                    if currentexon <= numExons:
                        if d['exonstart%s' % (str(currentexon))] <= startcodonstart and d['exonstop%s' % (str(currentexon))] >= startcodonstop:
                            if AFEstrand == '+':
                                startlist.append(startcodonstart) #start codon will be next 3 nt, useful for checking
                                break
                            elif AFEstrand == '-':
                                startlist.append(startcodonstop) #start codon will be next 3 nt, useful for checking
                                break
                        else:
                            currentexon +=1

        #Print number of unique start codons found for each AFE
        print len(list(set(startlist))) 

        if len(startlist) > 0:
            AFEswithStarts +=1
            if AFEstrand == '+':
                startpos = [AFEchrm, max(startlist)] #take most proximal start codon to give longest UTR
                isoformgff = [AFEchrm, 'AFE', '5UTR', str(AFEstart), str(startpos[1]), '.', AFEstrand, '.', AFEID]
                UTRgffAll.append(isoformgff)
                head, sep, tail = AFEID[3:].partition(';') #Remove 'ID=' and everything after semicolon
                for isoform in db.features_of_type('mRNA'):
                    if isoform.attributes['ID'] == head:
                        for exon in db.children(isoform, featuretype = 'exon'): #retrieve child exons
                            UTRgffAll.append([exon.chrom, 'AFE', exon.featuretype, str(exon.start), str(min(exon.stop, startpos[1])), '.', exon.strand, '.', 'ID=' + exon.attributes['ID'] + ';Parent=' + exon.attributes['Parent']]) #exon boundaries can't exceed UTR boundaries
                           

            elif AFEstrand == '-':
                startpos = [AFEchrm, min(startlist)]
                isoformgff = [AFEchrm, 'AFE', '5UTR', str(startpos[1]), str(AFEstop), '.', AFEstrand, '.', AFEID]
                UTRgffAll.append(isoformgff)
                head, sep, tail = AFEID[3:].partition(';') 
                for isoform in db.features_of_type('mRNA'):
                    if isoform.attributes['ID'] == head:
                        for exon in db.children(isoform, featuretype = 'exon'): 
                            UTRgffAll.append([exon.chrom, 'AFE', exon.featuretype, str(max(exon.start, startpos[1])), str(exon.stop), '.', exon.strand, '.', 'ID=' + exon.attributes['ID'] + ';Parent=' + exon.attributes['Parent']])

            
        #If no start codon found...
        elif len(startlist) == 0: 
            print 'WARNING: no start codon found for ' + AFEID 
        
        startlist = [] #reset startlist
        
    print 'Found start codons for %s of %s AFEs.' % (AFEswithStarts, len(AFEisoforms))
    return UTRgffAll
        
if __name__ == '__main__':
    outfh = open(sys.argv[3], 'w')
    for entry in get5UTRcoords(sys.argv[1], sys.argv[2]):
        outfh.write(('\t').join(entry) + '\n')

    outfh.close()

os.remove(db_fn)

    

    
