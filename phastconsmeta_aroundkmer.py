#Given a gff of regions you are interested in, a fasta file of the genome sequence, and a .bed file
#of phastcons values (sorted.phastcons.mm9.bed.gz), look for kmers within those regions and 
#get phastcons values from kmerstart - 100 to kmerend + 100

from Bio import SeqIO
import sys
import tabix
import pysam
from numpy import mean, std, sqrt
import argparse

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

def getkmerpos(gff, seq_dict, kmers):
    kmerpos = {} # {age : {location : [[chrm, kmerstart, kmerstop, strand]]}}
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
        ID = line[8][3:]
        region = (';').join([chrm, str(start), str(stop), strand])
        seq = getSequence(seq_dict, region)

        for i in range(len(seq) - k + 1):
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
                    kmerpos[age][location] = [[chrm, kmerstart, kmerstop, strand]]
                elif kmerpos.has_key(age):
                    if kmerpos[age].has_key(location) == False:
                        kmerpos[age][location] = [[chrm, kmerstart, kmerstop, strand]]
                    elif kmerpos[age].has_key(location):
                        kmerpos[age][location].append([chrm, kmerstart, kmerstop, strand])

    gfffh.close()
    return kmerpos

def getphastcons(kmerpos, phastconsbed):
    #kmerpos = {} # {age : {location : [[chrm, kmerstart, kmerstop, strand]]}}
    phastconsdict = {} # {age : {location : {binnumber : [phastconsvalue1, phastconsvalue2, ...]}}}
    #Bin1 is 100 bp upstream of kmerstart. For a 4mer, the kmer would be bins 101-104, and bins 105-205 would be 100 bp downstream of kmerstop.
    phastconstabix = pysam.Tabixfile(phastconsbed)
    for age in kmerpos:
        phastconsdict[age] = {}
        for location in kmerpos[age]:
            phastconsdict[age][location] = {}
            for kmer in kmerpos[age][location]:
                chrm = kmer[0]
                kmerstart = int(kmer[1])
                kmerstop = int(kmer[2])
                strand = kmer[3]
                phastconsscores = {} # {windowbin : score}
                if strand == '+':
                    windowstart = kmerstart - 100
                    windowend = kmerstop + 100
                    try:
                        for bed in phastconstabix.fetch(chrm, windowstart, windowend, parser = pysam.asBed()):
                            windowbin = str(int(bed.start) - windowstart)
                            phastconsscore = float(bed.name)
                            phastconsscores[windowbin] = phastconsscore
                    except ValueError:
                        print 'WARNING: problem with {0}:{1}-{2}:{3}.'.format(str(chrm), str(kmerstart), str(kmerstop), strand)

                elif strand == '-':
                    windowstart = kmerstart - 100
                    windowend = kmerstop + 100
                    try:
                        for bed in phastconstabix.fetch(chrm, windowstart, windowend, parser = pysam.asBed()):
                            windowbin = str(windowend - int(bed.start))
                            phastconsscore = float(bed.name)
                            phastconsscores[windowbin] = phastconsscore
                    except ValueError:
                        print 'WARNING: problem with {0}:{1}-{2}:{3}.'.format(str(chrm), str(kmerstart), str(kmerstop), strand)

                if len(phastconsscores) > 0: #if there were any bases in the UTR that had phastcons scores
                    for windowbin in phastconsscores:
                        if phastconsdict[age][location].has_key(windowbin) == False:
                            phastconsdict[age][location][windowbin] = [phastconsscores[windowbin]]
                        elif phastconsdict[age][location].has_key(windowbin):
                            phastconsdict[age][location][windowbin].append(phastconsscores[windowbin])

    return phastconsdict

def summarizephastconsdict(phastconsdict):
    #phastconsdict = {age : {location : {binnumber : [phastconsvalue1, phastconsvalue2, ...]}}}
    meanphastconsdict = {} # {age : {location : {binnumber : [meanphastcons, stdev]}}}
    for age in phastconsdict:
        meanphastconsdict[age] = {}
        for location in phastconsdict[age]:
            meanphastconsdict[age][location] = {}
            for binnumber in phastconsdict[age][location]:
                meanphastconsscore = round(mean(phastconsdict[age][location][binnumber]), 4)
                stdev = round(std(phastconsdict[age][location][binnumber]), 4)
                stderr = round(stdev / sqrt(len(phastconsdict[age][location][binnumber])), 4)
                meanphastconsdict[age][location][binnumber] = [meanphastconsscore, stdev, stderr]

    return meanphastconsdict

def smoothscores(meanphastconsdict):
	#meanphastconsdict = {age : {location : {binnumber : [meanphastcons, stdev]}}}
	smoothdict = {} # {age : {location : {binnumber : smoothedscore}}}
	smoothwindow = 2 #width of window to use when smoothing...this is the bin and +/- 2 nt
	#Bins range from 1 (100 bp upstream of kmer start) to 200 + k.  Kmer starts at 100

	for age in meanphastconsdict:
		if smoothdict.has_key(age) == False:
			smoothdict[age] = {}
		for location in meanphastconsdict[age]:
			if smoothdict[age].has_key(location) == False:
				smoothdict[age][location] = {}
			for binnumber in meanphastconsdict[age][location]:
				binnumber = int(binnumber)
				nearbyvalues = []
				if binnumber >=3 and binnumber <= 200: #cant look at the edges because there is no window
					for i in range(binnumber - smoothwindow, binnumber + smoothwindow + 1):
						if str(i) in meanphastconsdict[age][location]:
							nearbyvalues.append(meanphastconsdict[age][location][str(i)][0])
						else: #sometimes there is no value for a bin in a particular age/location
							continue
					smoothedvalue = mean(nearbyvalues)
					smoothdict[age][location][binnumber] = smoothedvalue

	return smoothdict


               
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
    parser.add_argument('--gff', type = str, help = 'Gff file of regions to consider.')
    parser.add_argument('--kmers', type = str, help = 'Comma separated list of kmers to look for.')
    parser.add_argument('--phastcons', type = str, help = 'Gzipped and tabix indexed bed of phastcons scores. sorted.phastcons.mm9.bed.gz works well. Requires sorted.phastcons.mm9.bed.gz.tbi in the same directory.')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    parser.add_argument('--smoothedoutfile', type = str, help = 'Optional. If you want the output smoothed, provide a filename here.')
    args = parser.parse_args()
    
    seq_dict = indexgenome(args.genomefasta)
    kmers = args.kmers.split(',')
    kmerpos = getkmerpos(args.gff, seq_dict, kmers)
    phastconsdict = getphastcons(kmerpos, args.phastcons)
    meanphastconsdict = summarizephastconsdict(phastconsdict)
    smoothdict = smoothscores(meanphastconsdict)

    def asint(s):
        try: 
            return int(s), ''
        
        except ValueError: 
            return sys.maxint, s
    
    outfh = open(args.outfile, 'w')
    outfh.write('Age' + '\t' + 'location' + '\t' + 'BNS' + '\t' + 'bin' + '\t' + 'meanphastcons' + '\t' + 'stdev' + '\t' + 'stderr' + '\t' + 'Protein' + '\t' + 'BNS' + '\n')
    for age in meanphastconsdict:
        for location in meanphastconsdict[age]:
            for binnumber in sorted(meanphastconsdict[age][location], key = asint):
                phastconsscore = meanphastconsdict[age][location][binnumber][0]
                stdev = meanphastconsdict[age][location][binnumber][1]
                stderr = meanphastconsdict[age][location][binnumber][2]
                outfh.write(age + '\t' + location + '\t' + 'foxunbound' + '\t' + str(binnumber) + '\t' + str(phastconsscore) + '\t' + str(stdev) + '\t' + str(stderr) + '\t' + 'Fox' + '\t' + 'unbound' + '\n')

    outfh.close()

    if args.smoothedoutfile:
    	outfh = open(args.smoothedoutfile, 'w')
    	outfh.write(('\t').join(['Age', 'location','Class', 'bin', 'smoothedphastcons','Protein','BNS']) + '\n')
    	for age in smoothdict:
    		for location in smoothdict[age]:
    			for binnumber in sorted(smoothdict[age][location], key = asint):
    				smoothedscore = smoothdict[age][location][binnumber]
    				outfh.write(('\t').join([age, location, 'msibound', str(binnumber), str(smoothedscore), 'Msi', 'bound']) + '\n')

    outfh.close()


