#This can work two ways, depending on the "mode".  In both modes, it calculates Rc for an oligo.  Rc is defined as R / Rexp, where Rexp is
#the median R for all oligos with that number of motifs.  In one mode, this returns the Rc for the oligo and the mean phastcons score for all
#bp within the oligo.  In another mode, this finds all kmer occurences within the oligos and returns phastcons scores +/- 10 bp around the 
#kmer as well as the R, Rexp, Rc values for the oligo that the kmer is contained within.


from Bio import SeqIO
import sys
import tabix
import pysam
import argparse
from numpy import mean, median

def indexgenome(genomefasta):
    sys.stderr.write('Indexing genome sequence...\n')
    seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
    sys.stderr.write('{0} chromosomes indexed.\n'.format(len(seq_dict)))

    return seq_dict

def getSequence(seq_dict, region):
    #region = chrm;start;stop;strand
    chrm = region.split(';')[0]
    start = int(region.split(';')[1])
    stop = int(region.split(';')[2])
    strand = region.split(';')[3]

    if strand == '+':
        seq = seq_dict[chrm].seq[start:stop].upper().transcribe()
    elif strand == '-':
        seq = seq_dict[chrm].seq[start:stop].upper().reverse_complement().transcribe()

    seq = str(seq)

    return seq

def getkmercount(gff, seq_dict, kmers):
    print 'Counting number of kmers in each sequence...'
    #Count the number of times one of the kmers occurs in every sequence
    kmercountdict = {} # {age : {location : {seqname : kmercount}}}
    k = len(kmers[0])
    gfffh = open(gff, 'r')
    for line in gfffh:
        kmercount = 0
        line = line.strip().split('\t')
        chrm = line[0]
        age = line[1]
        location = line[2]
        start = int(line[3])
        stop = int(line[4])
        strand = line[6]
        ID = line[8][3:]
        region = (';').join([chrm, str(start), str(stop), strand])
        seq = getSequence(seq_dict, region)

        #for i in range(len(seq) - k + 1):
            #kmer = seq[i:i+k]
            #if kmer in kmers:
                #kmercount +=1

     	i = 0
     	while i <= len(seq) - k + 1:
     		kmer = seq[i:i+k]
     		if kmer in kmers:
     			kmercount +=1
     			i +=k #no overlapping kmers
     		else:
     			i +=1

        if kmercountdict.has_key(age) == False:
            kmercountdict[age] = {}
        if kmercountdict[age].has_key(location) == False:
            kmercountdict[age][location] = {}
        kmercountdict[age][location][ID] = kmercount

    gfffh.close()

    return kmercountdict

def getRvalues(rvaluefile):
    print 'Retrieving R values...'
    #Retrieve R values for every sequence
    rdict = {} # {seqname (no ENSG) : rvalue}
    rfh = open(rvaluefile, 'r')
    for line in rfh:
        line = line.strip().split('\t')
        #skip header
        if line[0] == 'species':
            continue
        seqname = line[1]
        r = float(line[2])
        rdict[seqname] = r

    rfh.close()
    return rdict
    

def getRexpvalues(kmercountdict, rdict):
    print 'Calculating Rexp values...'
    #kmercountdict = {} # {age : {location : {seqname : kmercount}}}
    rexpdict = {} # {kmercount : rexp}
    kcount_r = {} # {kmercount : [rvalue1, rvalue2, ...]}
    for age in kmercountdict:
        for location in kmercountdict[age]:
            for seqname in kmercountdict[age][location]:
                kmercount = kmercountdict[age][location][seqname]
                ID = seqname.split('|')
                #Skip any berglund oligos
                if ID[0] == 'berglund':
                    continue
                newID = ('|').join([ID[0], ID[1], ID[2], ID[4], ID[5], ID[6], ID[7], ID[8]])
                rvalue = rdict[newID]
                if kcount_r.has_key(kmercount) == False:
                    kcount_r[kmercount] = [rvalue]
                elif kcount_r.has_key(kmercount) == True:
                    kcount_r[kmercount].append(rvalue)

    for kmercount in kcount_r:
        #averager = median(kcount_r[kmercount])
        #rexpdict[kmercount] = averager
        #Remove top and bottom 1% of R values
        rvalues = sorted(kcount_r[kmercount])
        bottom1percent = int(round(len(rvalues) * 0.01))
        top1percent = int(round(len(rvalues) * 0.99))
        averager = mean(kcount_r[kmercount][bottom1percent:top1percent])
        rexpdict[kmercount] = averager
    return rexpdict


def getphastcons(gff, phastconsbed):
    print 'Retrieving phastcons values...'
    #Get mean phastcons score for every oligo
    phastconsdict = {} # {seqname : meanphastcons}
    phastconstabix = pysam.Tabixfile(phastconsbed)
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        chrm = line[0]
        start = int(line[3])
        stop = int(line[4])
        strand = line[6]
        seqname = line[8][3:]
        scores = []
        try:
            for bed in phastconstabix.fetch(chrm, start, stop, parser = pysam.asBed()):
                phastconsscore = float(bed.name)
                scores.append(phastconsscore)
        except ValueError:
            print 'WARNING: problem with seq {0}, {1}:{2}-{3}.'.format(seqname, str(chrm), str(start), str(stop))
        if len(scores) > 0:
            meanphastcons = mean(scores)
            phastconsdict[seqname] = meanphastcons

    gfffh.close()
    return phastconsdict

def getmeanphastconsofregion(region, phastconsbed):
    #region = [chrm, start, stop, strand]
    chrm = region[0]
    start = int(region[1])
    stop = int(region[2])
    strand = region[3]
    scores = []
    phastconstabix = pysam.Tabixfile(phastconsbed)
    try:
        for bed in phastconstabix.fetch(chrm, start, stop, parser = pysam.asBed()):
            phastconsscore = float(bed.name)
            scores.append(phastconsscore)
    except ValueError:
        print 'WARNING: problem with seq {1}:{2}-{3}.'.format(str(chrm), str(start), str(stop))
    if len(scores) > 0:
        meanphastcons = mean(scores)
        return meanphastcons
    else:
        return None


def Rcandphastcons(rdict, rexpdict, phastconsdict, kmercountdict, gff):
    #This function outputs the Rc of an oligo and the mean phastcons score across its length
    #rdict = {} # {seqname (no ENSG) : rvalue}
    #rexpdict = {} # {kmercount : rexp}
    #phastconsdict = {} # {seqname : meanphastcons}
    #kmercountdict = {} # {age : {location : {seqname : kmercount}}}
    rcandphastconsdict = {} # {age : {location : {seqname : [rc, meanphastcons]}}}
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        age = line[1]
        location = line[2]
        seqname = line[8][3:]
        tempname = seqname.split('|')
        seqname_noensg = ('|').join([tempname[0], tempname[1], tempname[2], tempname[4], tempname[5], tempname[6], tempname[7], tempname[8]])
        
        #If it has a recorded R value and had at least one bp of phastcons score...
        if seqname_noensg in rdict and seqname in phastconsdict:
            rvalue = rdict[seqname_noensg]
            kmercount = kmercountdict[age][location][seqname]
            rexp = float(rexpdict[kmercount])
            rc = rvalue / rexp
            phastconsscore = phastconsdict[seqname]
            if rcandphastconsdict.has_key(age) == False:
                rcandphastconsdict[age] = {}
            if rcandphastconsdict[age].has_key(location) == False:
                rcandphastconsdict[age][location] = {}
            rcandphastconsdict[age][location][seqname] = [rc, phastconsscore]

    gfffh.close()
    return rcandphastconsdict

def getRcvalues(rdict, rexpdict, kmercountdict, gff):
    #Given a gff of oligos (just to get their names), its R value, the number of kmers it contains,
    #and the Rexp for that kmercount, calculate Rc
    #rdict = {} # {seqname (no ENSG) : rvalue}
    #rexpdict = {} # {kmercount : rexp}
    #kmercountdict = {} # {age : {location : {seqname : kmercount}}}
    rcdict = {} #{seqname : Rc}
    outfh = open('RcvsR.txt', 'w')
    outfh.write(('\t').join(['Age','Location','oligoname','kmercount','R','Rexp','Rc']) + '\n')
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        age = line[1]
        location = line[2]
        seqname = line[8][3:]
        tempname = seqname.split('|')
        seqname_noensg = ('|').join([tempname[0], tempname[1], tempname[2], tempname[4], tempname[5], tempname[6], tempname[7], tempname[8]])

        if seqname_noensg in rdict:
            rvalue = rdict[seqname_noensg]
            kmercount = kmercountdict[age][location][seqname]
            rexp = float(rexpdict[kmercount])
            rc = rvalue / rexp
            rcdict[seqname] = rc
            outfh.write(('\t').join([age, location, seqname, str(kmercount), str(rvalue), str(rexp), str(rc)]) + '\n')

    outfh.close()
    return rcdict


######STARTING HERE....LOOK AT Rc FOR AN OLIGO AND COMPARE IT TO PHASTCONS +/- 10 BP AROUND MOTIF IN OLIGO######
################################################################################################################

def getkmerpos(gff, seq_dict, kmers):
    #An oligo must have at least 1 instance of a kmer to end up in this dict
    kmerpos = {} # {age : {location : {oligoname : [[chrm, kmerstart, kmerstop, strand], [chrm, kmerstart, kmerstop, strand]]}}}
    k = len(kmers[0])
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        chrm = line[0]
        age = line[1]
        location = line[2]
        start = int(line[3])
        stop = int(line[4])
        strand = line[6]
        oligoname = line[8][3:]
        region = (';').join([chrm, str(start), str(stop), strand])
        seq = getSequence(seq_dict, region)

        #for i in range(len(seq) - k + 1):

        #NO overlapping kmers
        i = 0
        while i <= len(seq) - k + 1:
            kmer = seq[i:i+k]
            if kmer in kmers:
                if strand == '+':
                    kmerstart = start + i + 1
                    kmerstop = start + i + k
                elif strand == '-':
                    kmerstart = stop - i - k + 1
                    kmerstop = stop - i
                if kmerpos.has_key(age) == False:
                    kmerpos[age] = {}
                    kmerpos[age][location] = {}
                    kmerpos[age][location][oligoname] = [[chrm, kmerstart, kmerstop, strand]]
                elif kmerpos.has_key(age):
                    if kmerpos[age].has_key(location) == False:
                        kmerpos[age][location] = {}
                        kmerpos[age][location][oligoname] = [[chrm, kmerstart, kmerstop, strand]]
                    elif kmerpos[age].has_key(location):
                        if kmerpos[age][location].has_key(oligoname) == False:
                            kmerpos[age][location][oligoname] = [[chrm, kmerstart, kmerstop, strand]]
                        elif kmerpos[age][location].has_key(oligoname) == True:
                            kmerpos[age][location][oligoname].append([chrm, kmerstart, kmerstop, strand])

                i = i + k

            else:
                i += 1

    gfffh.close()
    return kmerpos

def kmerpostobed(kmerpos):
    outfh = open('Kmerpositions.bed', 'w')
    for age in kmerpos:
        for location in kmerpos[age]:
            for oligoname in kmerpos[age][location]:
                for kmerinstance in kmerpos[age][location][oligoname]:
                    chrm = kmerinstance[0]
                    kmerstart = str(kmerinstance[1] - 1) #bed is 0-based
                    kmerstop = str(kmerinstance[2]) #is this right?  [start : stop)
                    strand = kmerinstance[3]
                    outfh.write(('\t').join([chrm, kmerstart, kmerstop, oligoname, '1000', strand]) + '\n')

    outfh.close()


def getphastconsaroundkmer(kmerpos, rdict, rexpdict, phastconsbed):
    #Given a set of kmer positions and Rexp values for every oligo, calculate mean phastcons score
    #+/- 10 bp around kmer
    #kmerpos = {age : {location : {oligoname : [[chrm, kmerstart, kmerstop, strand], [chrm, kmerstart, kmerstop, strand]]}}}
    #rdict = {} # {seqname (no ENSG) : rvalue}
    #rexpdict = {} # {kmercount : rexp}
    rckmerdict = {} #{age : {location : {oligoname : [[chrm, kmerstart, kmerstop, strand, r, rexp, rc, phastcons], [chrm, kmerstart, kmerstop, strand, r, rexp, rc, phastcons]]}}}


    for age in kmerpos:
        if rckmerdict.has_key(age) == False:
            rckmerdict[age] = {}
        for location in kmerpos[age]:
            print 'Finding phastcons scores around kmers in {0} {1} oligos.'.format(age, location)
            if rckmerdict[age].has_key(location) == False:
                rckmerdict[age][location] = {}
            for oligoname in kmerpos[age][location]:
                if rckmerdict[age][location].has_key(oligoname) == False:
                    #make an empty list
                    rckmerdict[age][location][oligoname] = []
                for i in range(len(kmerpos[age][location][oligoname])): #Here, every i is a kmer occurence
                    chrm = kmerpos[age][location][oligoname][i][0]
                    kmerstart = kmerpos[age][location][oligoname][i][1]
                    kmerstop = kmerpos[age][location][oligoname][i][2]
                    strand = kmerpos[age][location][oligoname][i][3]
                    tempname = oligoname.split('|')
                    if tempname[0] == 'berglund': #skip any berglund oligos
                        continue
                    oligoname_noensg = ('|').join([tempname[0], tempname[1], tempname[2], tempname[4], tempname[5], tempname[6], tempname[7], tempname[8]])
                    rvalue = rdict[oligoname_noensg]
                    numberofkmersinoligo = len(kmerpos[age][location][oligoname])
                    rexp = rexpdict[numberofkmersinoligo] 
                    rc = rvalue / rexp
                    regionaroundkmer = [chrm, kmerstart-10, kmerstop+10, strand] #we will define the region of interest around the kmer as +/- 10 bp
                    meanphastcons = getmeanphastconsofregion(regionaroundkmer, phastconsbed)
                    if meanphastcons:
                        #if there were no bp with phastcons score in the region, meanphastcons == None
                        rckmerdict[age][location][oligoname].append([chrm, kmerstart, kmerstop, strand, numberofkmersinoligo, rvalue, rexp, rc, meanphastcons])

    return rckmerdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
    parser.add_argument('--kmers', type = str, help = 'Comma separated list of kmers to look for.')
    parser.add_argument('--gff', type = str, help = 'Gff file of oligo regions to consider.')
    parser.add_argument('--rvalues', type = str, help = 'R value file: tab separated species, oligoname, rvalue, motif_num. Only oligoname and rvalue are used.')
    parser.add_argument('--phastconsscores', type = str, help = 'Gzipped and tabix indexed bed of phastcons scores. sorted.phastcons.mm9.bed.gz works well. Requires sorted.phastcons.mm9.bed.gz.tbi in the same directory.')
    parser.add_argument('--mode', type = str, choices = ['scorewholeoligo', 'scorearoundkmer', 'getRcvalues', 'kmerposbed'], help = 'Do you want phastconsscores of the whole oligo or phastconsscores just around kmers or Rc values of oligos or a bed file of kmer positions? If Rc values, output will be in RvsRc.txt. If bed of kmer positions, output is Kmerpositions.bed')
    parser.add_argument('--output', type = str, help = 'Output file.')
    args = parser.parse_args()

    if args.mode == 'scorewholeoligo':

        seq_dict = indexgenome(args.genomefasta)
        kmers = args.kmers.split(',')
        kmercountdict = getkmercount(args.gff, seq_dict, kmers)
        rdict = getRvalues(args.rvalues)
        rexpdict = getRexpvalues(kmercountdict, rdict)
        print rexpdict
        phastconsdict = getphastcons(args.gff, args.phastconsscores)
        rcandphastconsdict = Rcandphastcons(rdict, rexpdict, phastconsdict, kmercountdict, args.gff)
        outfh = open(args.output, 'w')
        outfh.write(('\t').join(['Age','Location','sequence','Rc','MeanPhastcons']) + '\n')
        for age in rcandphastconsdict:
            for location in rcandphastconsdict[age]:
                for seqname in rcandphastconsdict[age][location]:
                    rc = rcandphastconsdict[age][location][seqname][0]
                    meanphastcons = rcandphastconsdict[age][location][seqname][1]
                    outfh.write(age + '\t' + location + '\t' + seqname + '\t' + str(round(rc, 4)) + '\t' + str(round(meanphastcons, 4)) + '\n')

        outfh.close()

    elif args.mode == 'scorearoundkmer':
        seq_dict = indexgenome(args.genomefasta)
        kmers = args.kmers.split(',')
        kmerpos = getkmerpos(args.gff, seq_dict, kmers)
        kmercountdict = getkmercount(args.gff, seq_dict, kmers)
        rdict = getRvalues(args.rvalues)
        rexpdict = getRexpvalues(kmercountdict, rdict)
        rckmerdict = getphastconsaroundkmer(kmerpos, rdict, rexpdict, args.phastconsscores)
        outfh = open(args.output, 'w')
        #rckmerdict = {age : {location : {oligoname : [[chrm, kmerstart, kmerstop, strand, r, rexp, rc, phastcons], [chrm, kmerstart, kmerstop, strand, r, rexp, rc, phastcons]]}}}
        outfh.write(('\t').join(['Age','Location','Oligoname','chrm','kmerstart','kmerstop','strand', 'kmers_in_oligo','r','rexp','rc','mean_phastcons_around_kmer']) + '\n')
        for age in rckmerdict:
            for location in rckmerdict[age]:
                for oligoname in rckmerdict[age][location]:
                    for kmerinstance in rckmerdict[age][location][oligoname]:
                        chrm = kmerinstance[0]
                        kmerstart = str(kmerinstance[1])
                        kmerstop = str(kmerinstance[2])
                        strand = kmerinstance[3]
                        kmers_in_oligo = str(kmerinstance[4])
                        r = str(kmerinstance[5])
                        rexp = str(kmerinstance[6])
                        rc = str(kmerinstance[7])
                        meanphastcons = str(kmerinstance[8])
                        outfh.write(('\t').join([age, location, oligoname, chrm, kmerstart, kmerstop, strand, kmers_in_oligo, r, rexp, rc, meanphastcons]) + '\n')

        outfh.close()

    elif args.mode == 'getRcvalues':
        seq_dict = indexgenome(args.genomefasta)
        kmers = args.kmers.split(',')
        kmercountdict = getkmercount(args.gff, seq_dict, kmers)
        rdict = getRvalues(args.rvalues)
        rexpdict = getRexpvalues(kmercountdict, rdict)
        getRcvalues(rdict, rexpdict, kmercountdict, args.gff)

    elif args.mode == 'kmerposbed':
        seq_dict = indexgenome(args.genomefasta)
        kmers = args.kmers.split(',')
        kmerpos = getkmerpos(args.gff, seq_dict, kmers)
        kmerpostobed(kmerpos)
