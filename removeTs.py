#python 3
#usage python removeTs.py <infile.gz> <outfile.gz>

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
import gzip

#Remove Ts from the beginning of sequencing reads, like if you have 3' end data

#make output file
with open(sys.argv[2], 'w') as outfh:
	pass

counter = 0
with open(sys.argv[1], 'r') as infh, open(sys.argv[2], 'a') as outfh:
	for title, seq, qual in FastqGeneralIterator(infh):
		counter +=1
		if counter % 1000000 == 0:
			print('Read {0}...'.format(counter))
		if not seq.startswith('TT'):
			continue
		else:
			if seq.startswith('TTTTTTTTTTTT'):
				seq = seq[12:]
				qual = qual[12:]

			elif seq.startswith('TTTTTTTTTTT'):
				seq = seq[11:]
				qual = qual[11:]

			elif seq.startswith('TTTTTTTTTT'):
				seq = seq[10:]
				qual = qual[10:]					

			elif seq.startswith('TTTTTTTTT'):
				seq = seq[9:]
				qual = qual[9:]

			elif seq.startswith('TTTTTTTT'):
				seq = seq[8:]
				qual = qual[8:]

			elif seq.startswith('TTTTTTT'):
				seq = seq[7:]
				qual = qual[7:]

			elif seq.startswith('TTTTTT'):
				seq = seq[6:]
				qual = qual[6:]

			elif seq.startswith('TTTTT'):
				seq = seq[5:]
				qual = qual[5:]

			elif seq.startswith('TTTT'):
				seq = seq[4:]
				qual = qual[4:]

			elif seq.startswith('TTT'):
				seq = seq[3:]
				qual = qual[3:]

			elif seq.startswith('TT'):
				seq = seq[2:]
				qual = qual[2:]

			outfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq, qual))




