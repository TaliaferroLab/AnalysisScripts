import subprocess
import argparse
import sys
import os

def makeindex(transcripts):
	#Transcripts is fasta file of transcripts to be quantified.
	#Can be gzipped.

	command = ['kallisto', 'index', '-k', '31', '-i', 'transcripts.idx', transcripts]
	subprocess.call(command)

def runkallisto(outputname, threads, reads1, reads2):
	#From paired end sequencing data

	#Paired end
	if reads2:
		command = ['kallisto', 'quant', '-i', 'transcripts.idx', '-o', outputname, '-b', '100', '--bias', '-t', threads, reads1, reads2]

	#Single end, setup here for ribosome footprint data
	elif not reads2:
		s = '10'
		if 'rpf' in outputname:
			l = '30'
			command = ['kallisto', 'quant', '-i', 'transcripts19.idx', '-o', outputname, '-b', '100', '--bias', '--single', '-l', l, '-s', s, '-t', threads, reads1]
		elif 'tot' in outputname:
			l = '40'
			command = ['kallisto', 'quant', '-i', 'transcripts25.idx', '-o', outputname, '-b', '100', '--bias', '--single', '-l', l, '-s', s, '-t', threads, reads1]

	subprocess.call(command)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--transcripts', type = str, help = 'Fasta file of transcripts to quantify.')
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
	makeindex(args.transcripts)
	print 'Done indexing!'

	if args.reversereads:
		for i in range(len(labels)):
			label = labels[i]
			freads = forwardreads[i]
			rreads = reversereads[i]
			print 'Quantifying sample {0}...'.format(label)
			runkallisto(label, threads, freads, rreads)

	elif not args.reversereads:
		for i in range(len(labels)):
			label = labels[i]
			freads = forwardreads[i]
			rreads = None
			print 'Quantifying sample {0}...'.format(label)
			runkallisto(label, threads, freads, rreads)

os.remove('transcripts.idx')

