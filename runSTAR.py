import os
import subprocess
import argparse
import sys

#python3

def runSTAR_paired(read1, read2, samplename, outputdir, genomedir):
	print('Running STAR for {0}...'.format(samplename))

	resultsdir = os.path.join(os.path.abspath(outputdir), samplename)
	os.mkdir(resultsdir)
	prefix = resultsdir + '/' + samplename
	command = ['STAR', '--runMode', 'alignReads', '--runThreadN', '8', '--genomeLoad', 'NoSharedMemory', '--genomeDir', genomedir, '--readFilesIn', read1, read2, '--readFilesCommand', 'zcat',
	'--outFileNamePrefix', prefix, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', 'intronMotif']

	subprocess.call(command)

	print('Finished STAR for {0}!'.format(samplename))

def runSTAR_single(read1, samplename, outputdir, genomedir):
	print('Running STAR for {0}...'.format(samplename))

	resultsdir = os.path.join(os.path.abspath(outputdir), samplename)
	os.mkdir(resultsdir)
	prefix = resultsdir + '/' + samplename
	command = ['STAR', '--runMode', 'alignReads', '--runThreadN', '8', '--genomeLoad', 'NoSharedMemory', '--genomeDir', genomedir, '--readFilesIn', read1, '--readFilesCommand', 'zcat',
	'--outFileNamePrefix', prefix, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', 'intronMotif']

	subprocess.call(command)

	print('Finished STAR for {0}!'.format(samplename))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--readtype', choices = ['single', 'paired'], type = str, help = 'Single or paired end reads?')
	parser.add_argument('--forwardreads', type = str, help = 'Comma separated list of forward reads.')
	parser.add_argument('--reversereads', type = str, help = 'Comma separated list of reverse reads. Not needed if using single end reads.')
	parser.add_argument('--samplenames', type = str, help = 'Comma separated list of sample names.')
	parser.add_argument('--outputdir', type = str, help = 'Output directory where all the sample outputs will go.')
	parser.add_argument('--genomedir', type = str, help = 'Path to STAR index.')
	args = parser.parse_args()

	if args.readtype == 'paired':
		forwardreads = args.forwardreads.split(',')
		reversereads = args.reversereads.split(',')
		samplenames = args.samplenames.split(',')

		if len(forwardreads) != len(reversereads) != len(samplenames):
			print('ERROR: lists of read files and sample names must be same length!')
			sys.exit()

		for i in range(len(forwardreads)):
			read1 = forwardreads[i]
			read2 = reversereads[i]
			samplename = samplenames[i]
			runSTAR_paired(read1, read2, samplename, args.outputdir, args.genomedir)

	elif args.readtype == 'single':
		forwardreads = args.forwardreads.split(',')
		samplenames = args.samplenames.split(',')

		if len(forwardreads) != len(samplenames):
			print('ERROR: lists of read files and sample names must be same length!')
			sys.exit()

		for i in range(len(forwardreads)):
			read1 = forwardreads[i]
			samplename = samplenames[i]
			runSTAR_single(read1, samplename, args.outputdir, args.genomedir)
