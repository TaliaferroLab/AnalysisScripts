from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess
import argparse
import gzip
import sys
import os

#python2


def makeindex(transcripts, k):
	#Transcripts is fasta file of transcripts to be quantified.
	#Can be gzipped.

	k = str(k)

	command = ['salmon', 'index', '-t', transcripts, '-i', 'transcripts.idx', '--type', 'quasi', '-k', k]
	subprocess.call(command)

def runsalmon(outputname, threads, reads1, reads2):

	#paired end
	if reads2:
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '--validateMappings', '-1', reads1, '-2', reads2, '-o', outputname, '--index', 'transcripts.idx']

	#Single end, setup here for ribosome footprint data
	elif not reads2:
		fldMean = '35' #fragment length distribution mean
		fldSD = '10' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', outputname, '--index', 'transcripts.idx']

	subprocess.call(command)

def fastqlengthfilter(allowedlengths, infile):
	#Maybe you only want to allow reads of certain lengths
	#This is currently only set up to handle single end reads.  It's meant for something like ribosome footprinting.

	allowedlengths = [int(l) for l in allowedlengths]

	counter = 0
	passingreads = 0
	outfilename = 'temp.fastq'
	with gzip.open(infile, 'rb') as infh, open(outfilename, 'w') as outfh:
		try:
			for title, seq, qual in FastqGeneralIterator(infh):
				counter +=1
				if counter % 100000000 == 0:
					print 'On read {0} of {1}.'.format(counter, infile)
				if len(seq) in allowedlengths:
					passingreads +=1
					outfh.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

		except ValueError:
			print 'Title and second title line don\'t match for read {0}.'.format(title)

	print '{0} of {1} reads ({2} %) pass length filters.'.format(passingreads, counter, round((passingreads / float(counter)) * 100, 4))

def fastqtrimmer(threeprimetrim, forreads, revreads):
	#Maybe you want to trim the reads from the 3' end before giving them to salmon.
	#Trim <threeprimetrim> nt from the 3' end of the read
	threeprimetrim = int(threeprimetrim)

	counter = 0
	foutfilename = 'tempf.fastq'
	routfilename = 'tempr.fastq'
	with gzip.open(forreads, 'rb') as forinfh, gzip.open(revreads, 'rb') as revinfh, open(foutfilename, 'w') as foroutfh, open(routfilename, 'w') as revoutfh:
		try:
			for title, seq, qual in FastqGeneralIterator(forinfh):
				counter +=1
				if counter % 1000000 == 0:
					print 'On read {0} of {1}.'.format(counter, forreads)
				foroutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))
		except ValueError:
			pass

		try:
			for title, seq, qual in FastqGeneralIterator(revinfh):
				counter +=1
				if counter % 1000000 == 0:
					print 'On read {0} of {1}.'.format(counter, forreads)
				revoutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))

		except ValueError:
			pass

	print 'Done trimming {0} and {1}.'.format(forreads, revreads)





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--transcripts', type = str, help = 'Fasta file of transcripts to quantify.')
	parser.add_argument('-k', type = int, help = 'Value of k to use when making the index.')
	parser.add_argument('--allowedlengths', type = str, help = 'Optional. Comma separated list of allowed read lengths.')
	parser.add_argument('--threeprimetrim', type = int, help = 'Optional. Number of bases to trim from the 3\' end of reads.')
	parser.add_argument('--labels', type = str, help = 'Comma separated list of labels for experiments.')
	parser.add_argument('--threads', type = str)
	parser.add_argument('--forwardreads', type = str, help = 'Comma separated list of forward fastq files in same order as labels.')
	parser.add_argument('--reversereads', type = str, help = 'Comma separated list of reverse fastq files in same order as labels. Only needed with paired end reads.')
	args = parser.parse_args()

	labels = args.labels.split(',')
	forwardreads = args.forwardreads.split(',')
	if args.reversereads:
		reversereads = args.reversereads.split(',')
	threads = args.threads

	if args.reversereads:
		if len(labels) != len(forwardreads) != len(reversereads):
			print 'Lengths of labels and/or reads are not equal!'
			sys.exit()

	elif not args.reversereads:
		if len(labels) != len(forwardreads):
			print 'Lengths of labels and reads are not equal!'
			sys.exit()

	print 'Indexing transcripts...'
	makeindex(args.transcripts, args.k)
	print 'Done indexing!'

	if args.reversereads:
		if not args.threeprimetrim:
			for i in range(len(labels)):
				label = labels[i]
				freads = forwardreads[i]
				rreads = reversereads[i]
				print 'Quantifying sample {0}...'.format(label)
				runsalmon(label, threads, freads, rreads)

		elif args.threeprimetrim:
			for i in range(len(labels)):
				label = labels[i]
				freads = forwardreads[i]
				rreads = reversereads[i]
				print 'Trimming sample {0}...'.format(label)
				fastqtrimmer(args.threeprimetrim, freads, rreads)
				print 'Quantifying sample {0}...'.format(label)
				runsalmon(label, threads, 'tempf.fastq', 'tempr.fastq')
			os.remove('tempf.fastq')
			os.remove('tempr.fastq')

	elif not args.reversereads:
		if not args.allowedlengths:
			for i in range(len(labels)):
				label = labels[i]
				freads = forwardreads[i]
				rreads = None
				print 'Quantifying sample {0}...'.format(label)
				runsalmon(label, threads, freads, rreads)

		elif args.allowedlengths:
			allowedlengths = args.allowedlengths.split(',')
			for i in range(len(labels)):
				label = labels[i]
				freads = forwardreads[i]
				rreads = None
				print 'Trimming sample {0}...'.format(label)
				fastqlengthfilter(allowedlengths, freads)
				runsalmon(label, threads, 'temp.fastq', rreads)
				os.remove('temp.fastq')

