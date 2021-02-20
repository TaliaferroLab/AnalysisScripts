#For each event in rMATS output files, get the 5' splice site sequences.
#This really only makes sense for A5SS (there will be two sites per event),
#SE (there will be two sites per event), and RI (there will be one site per event)

from Bio import SeqIO
import os
import sys
import gzip

def makeseqdict(genomefasta):
    print('Indexing genome...')
    with gzip.open(genomefasta, 'rt') as infh:
        seq_dict = SeqIO.to_dict(SeqIO.parse(infh, 'fasta'))
    print('Done!')

    return seq_dict

def get5ssforA5ss(A5ssrmats, seq_dict, outfile):
    with open(A5ssrmats, 'r') as infh, open(outfile, 'w') as outfh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'ID':
                continue
            gene = line[1][1:-1] #remove quotes
            chrm = line[3]
            strand = line[4]
            longexonstart = int(line[5])
            longexonend = int(line[6])
            shortexonstart = int(line[7])
            shortexonend = int(line[8])
            ID = ':'.join([chrm, strand, str(longexonstart), str(longexonend), str(shortexonstart), str(shortexonend)])
            #Shorter isoform (one that uses upstream 5' ss is exclusion isoform)
            #5ss is 3 nt of exon and 6 nt of intron
            if strand == '+':
                exclusionssstart = shortexonend - 3
                exclusionssend = shortexonend + 6
                exclusionss = seq_dict[chrm].seq[exclusionssstart : exclusionssend].upper()
                inclusionssstart = longexonend - 3
                inclusionssend = longexonend + 6
                inclusionss = seq_dict[chrm].seq[inclusionssstart : inclusionssend].upper()
            elif strand == '-':
                exclusionssstart = shortexonstart - 6
                exclusionssend = shortexonstart + 3
                exclusionss = seq_dict[chrm].seq[exclusionssstart : exclusionssend].reverse_complement().upper()
                inclusionssstart = longexonstart - 6
                inclusionssend = longexonstart + 3
                inclusionss = seq_dict[chrm].seq[inclusionssstart : inclusionssend].reverse_complement().upper()
            outfh.write('>' + ID + ':exclusion' + '\n' + str(exclusionss) + '\n' + '>' + ID + ':inclusion' + '\n' + str(inclusionss) + '\n')

def get5ssforSE(SErmats, seq_dict, outfile):
    with open(SErmats, 'r') as infh, open(outfile, 'w') as outfh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'ID':
                continue
            gene = line[1][1:-1]
            chrm = line[3]
            strand = line[4]
            skippedES = int(line[5])
            skippedEE = int(line[6])
            upstreamES = int(line[7])
            upstreamEE = int(line[8])
            downstreamES = int(line[9])
            downstreamEE = int(line[10])
            ID = ':'.join([chrm, strand, str(skippedES), str(skippedEE), str(upstreamES), str(upstreamEE), str(downstreamES), str(downstreamEE)])
            if strand == '+':
                exclusionssstart = upstreamEE - 3
                exclusionssend = upstreamEE + 6
                exclusionss = seq_dict[chrm].seq[exclusionssstart : exclusionssend].upper()
                inclusionssstart = skippedEE - 3
                inclusionssend = skippedEE + 6
                inclusionss = seq_dict[chrm].seq[inclusionssstart : inclusionssend].upper()
            elif strand == '-':
                exclusionssstart = upstreamES - 6
                exclusionssend = upstreamES + 3
                exclusionss = seq_dict[chrm].seq[exclusionssstart : exclusionssend].reverse_complement().upper()
                inclusionssstart = skippedES - 6 
                inclusionssend = skippedES + 3
                inclusionss = seq_dict[chrm].seq[inclusionssstart : inclusionssend].reverse_complement().upper()
            outfh.write('>' + ID + ':exclusion' + '\n' + str(exclusionss) + '\n' + '>' + ID + ':inclusion' + '\n' + str(inclusionss) + '\n')

def get5ssforRI(RIrmats, seq_dict, outfile):
    with open(RIrmats, 'r') as infh, open(outfile, 'w') as outfh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'ID':
                continue
            gene = line[1][1:-1]
            chrm = line[3]
            strand = line[4]
            ristart = int(line[5])
            riend = int(line[6])
            upstreamES = int(line[7])
            upstreamEE = int(line[8])
            downstreamES = int(line[9])
            downstreamEE = int(line[10])
            ID = ':'.join([chrm, strand, str(ristart), str(riend), str(upstreamES), str(upstreamEE), str(downstreamES), str(downstreamEE)])
            if strand == '+':
                ssstart = upstreamEE - 3
                ssend = upstreamEE + 6
                ss = seq_dict[chrm].seq[ssstart : ssend].upper()
            elif strand == '-':
                ssstart = upstreamES - 6
                ssend = upstreamES + 3
                ss = seq_dict[chrm].seq[ssstart : ssend].reverse_complement().upper()
            outfh.write('>' + ID + '\n' + str(ss) + '\n')

#OK so after supplying these fastas to MaxEnt, you get back an output that needs to be parsed.

#This can parse both a5 output and SE output
def parsemaxenta5(a5maxentout, outfile):
    with open(a5maxentout, 'r') as infh, open(outfile, 'w') as outfh:
        outfh.write(('\t').join(['ID', 'exclusionss', 'exclusionscore', 'inclusionss', 'inclusionscore']) + '\n')
        for line in infh:
            line = line.strip().split()
            if line[0].startswith('>'):
                ID = line[0][1:-10]
                #exclusion always comes first
                if line[0].endswith('exclusion'):
                    site = 'exclusion'
                elif line[0].endswith('inclusion'):
                    site = 'inclusion'
            else:
                if site == 'exclusion':
                    #this is an exclusion score line
                    exclusionss = line[0]
                    exclusionscore = line[2]
                elif site == 'inclusion':
                    inclusionss = line[0]
                    inclusionscore = line[2]
                    outfh.write(('\t').join([ID, exclusionss, exclusionscore, inclusionss, inclusionscore]) + '\n')
                    exclusionss = None
                    exclusionscore = None
                    inclusionss = None
                    inclusionscore = None


def parsemaxentri(rimaxentout, outfile):
    with open(rimaxentout, 'r') as infh, open(outfile, 'w') as outfh:
        outfh.write(('\t').join(['ID', 'ss', 'ssscore']) + '\n')
        for line in infh:
            line = line.strip().split()
            if line[0].startswith('>'):
                ID = line[0][1:]
            else:
                ss = line[0]
                ssscore = line[2]
                outfh.write(('\t').join([ID, ss, ssscore]) + '\n')

#usage: python rMATS_5ss.py <gzipped genomefasta> <rMATSoutput> <outfile>

#seq_dict = makeseqdict(sys.argv[1])
#get5ssforA5ss(sys.argv[2], seq_dict, sys.argv[3])
#get5ssforSE(sys.argv[2], seq_dict, sys.argv[3])
#get5ssforRI(sys.argv[2], seq_dict, sys.argv[3])

#For maxent parsing
#parsemaxenta5(sys.argv[1], sys.argv[2])
parsemaxentri(sys.argv[1], sys.argv[2])