#Given a gff of regions you are interested in, a fasta file of the genome sequence, and a .bed file
#of phastcons values (sorted.phascons.mm9.bed.gz), look for kmers within those regions and get 
#phastconds values from kmerstart -25 to kmerend + 25.  Then take the mean of the region as the score
#for that oligo

from Bio import SeqIO
import sys
import tabix
import pysam
from numpy import mean, std, sqrt
import argparse
from phastconsmeta_aroundkmer import indexgenome, getSequence, getkmerpos
import os


def getphastcons(kmerpos, phastconsbed, outfile, protein, RBNSstate):
	#kmerpos = {} # {age : {location : [[chrm, kmerstart, kmerstop, strand]]}}
    phastconsdict = {} # {age : {location : [meanphastcons of oligo1 around motif, meanphastcons of oligo2 around motif]}}
    phastconstabix = pysam.Tabixfile(phastconsbed)
    for age in kmerpos:
    	phastconsdict[age] = {}
    	for location in kmerpos[age]:
    		phastconsdict[age][location] = []
    		for kmer in kmerpos[age][location]:
    			chrm = kmer[0]
    			kmerstart = int(kmer[1])
    			kmerstop = int(kmer[2])
    			strand = kmer[3]
    			phastconsscores = []
    			windowstart = kmerstart - 25
    			windowend = kmerstop + 25
    			try:
    				for bed in phastconstabix.fetch(chrm, windowstart, windowend, parser = pysam.asBed()):
    					phastconsscore = float(bed.name)
    					phastconsscores.append(phastconsscore)
    			except ValueError:
    				print 'WARNING: problem with {0}:{1}-{2}:{3}.'.format(str(chrm), str(kmerstart), str(kmerstop), strand)

    			if len(phastconsscores) > 0: #if there were any bases in the region that had phastcons scores
    				meanphastcons = mean(phastconsscores)
    				phastconsdict[age][location].append(meanphastcons)

    if not os.path.isfile(outfile):
    	with open(outfile, 'w') as f:
    		f.write(('\t').join(['age', 'location', 'protein', 'RBNSstate', 'meanphastcons']) + '\n')

    for age in phastconsdict:
    	for location in phastconsdict[age]:
    		for score in phastconsdict[age][location]:
    			with open(outfile, 'a') as f:
    				f.write(('\t').join([age, location, protein, RBNSstate, str(score)]) + '\n')




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--gff', type = str, help = 'Gff file of regions to consider.')
	parser.add_argument('--kmers', type = str, help = 'Comma separated list of kmers to look for.')
	parser.add_argument('--phastconsbed', type = str, help = 'Gzipped and tabix indexed bed of phastcons scores. sorted.phastcons.mm9.bed.gz works well. Requires sorted.phastcons.mm9.bed.gz.tbi in the same directory.')
	parser.add_argument('--protein', type = str, help = 'Name of protein that corresponds to the kmers.')
	parser.add_argument('--RBNSstate', type = str, help = 'RBNS state of the oligos encoded by the gff.  Usually \'bound\' or \'unbound\'.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')

	args = parser.parse_args()
	seq_dict = indexgenome(args.genomefasta)
	kmerpos = getkmerpos(args.gff, seq_dict, args.kmers)
	getphastcons(kmerpos, args.phastconsbed, args.outfile, args.protein, args.RBNSstate)
