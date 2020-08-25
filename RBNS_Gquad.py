#Usage: python RBNS_Gquad.py <fastqfile> <outfile>

from Bio import SeqIO
import gzip
import re
import itertools
import sys

def getdensity(motif, fastq):
	readcount = 0
	motifcount = 0
	fh = gzip.open(fastq, 'r')
	for record in SeqIO.parse(fh, 'fastq'):
		readcount +=1
		if readcount % 1000000 == 0:
			print 'Analyzing read {0}...'.format(readcount)
		seq = str(record.seq)
		motifmatches = re.findall(motif, seq)
		motifcount += len(motifmatches)

	return readcount, motifcount

def iteratemotifs(fastq, outfile):
	with open(outfile, 'a') as f:
		f.write(('\t').join(['fourmer', 'readcount', 'motifcount']) + '\n')
	fourmercount = 0
	bases = ['A', 'T', 'G', 'C']
	all4mers = [''.join(x) for x in itertools.product(bases, repeat = 4)]
	for fourmer in all4mers:
		fourmercount +=1
		if fourmercount % 10 == 0:
			print '4mer {0} of 256...'.format(fourmercount)
		motif = r'(?=({0}(.{{0,5}}){0}(.{{0,5}}){0}(.{{0,5}}){0}))'.format(fourmer)
		readcount, motifcount = getdensity(motif, fastq)
		with open(outfile, 'a') as f:
			f.write(('\t').join([fourmer, str(readcount), str(motifcount)]) + '\n')


iteratemotifs(sys.argv[1], sys.argv[2])