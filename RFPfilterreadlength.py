from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
import os
import gzip
import argparse

def fastqtrimmer(allowedlengths, infile):
	counter = 0
	outfilename = infile.split('.')[0] + '.' + infile.split('.')[1] + '.lengthfiltered.fastq.gz'
	with gzip.open(infile, 'rb') as infh, gzip.open(outfilename, 'wb') as outfh:
		try:
			for title, seq, qual in FastqGeneralIterator(infh):
				counter +=1
				if counter % 1000000 == 0:
					print 'On read {0} of {1}.'.format(counter, infile)
				if len(seq) in allowedlengths:
					outfh.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

		except ValueError:
			print 'Title and second title line don\'t match for read {0}.'.format(title)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fastqfiles', type = str, help = 'Comma separated list of fastq files to trim.')
	parser.add_argument('--allowedlengths', type = str, help = 'Comma separated list of allowed read lengths.')
	args = parser.parse_args()

	allowedlengths = args.allowedlengths.split(',')
	allowedlengths = [int(allowedlength) for allowedlength in allowedlengths]
	fastqfiles = args.fastqfiles.split(',')
	for fastqfile in fastqfiles:
		fastqtrimmer(allowedlengths, fastqfile)



