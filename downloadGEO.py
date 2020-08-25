import subprocess
import os
import sys

def downloadSRA(srr):
	#Downloads files to /vol3/home/taliaferro/ncbi/public/sra
	command = ['prefetch', '--max-size', '100G', '-v', srr]
	subprocess.call(command)

def sra2fastq_paired(sra, outdir):
	command = ['fasterq-dump', '-p', '--outdir', outdir, sra]
	subprocess.call(command)

def sra2fastq_single(sra, outdir):
	command = ['fasterq-dump', '-p', '--outdir', outdir, sra]
	subprocess.call(command)


def dllist(dlfile, end):
	srrfiles = []
	with open(dlfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			srrfiles.append(line[0])

	for idx, srr in enumerate(srrfiles):
		print 'Processing {0}, file {1} of {2}...'.format(srr, idx + 1, len(srrfiles))

		downloadSRA(srr)
		sra = os.path.join('/beevol/home/taliaferro/ncbi/public/sra', srr + '.sra')
		if end == 'single':
			sra2fastq_single(sra, '/beevol/home/taliaferro/data/TXpsi/Farris/RawReads')
			#Gzip files
			command = ['gzip', srr + '.sra.fastq']
			subprocess.call(command)
		elif end == 'paired':
			sra2fastq_paired(sra, '/beevol/home/taliaferro/data/TXpsi/Farris/RawReads')
			#Gzip files
			command = ['gzip', srr + '_1.fastq']
			subprocess.call(command)
			command = ['gzip', srr + '_2.fastq']
			subprocess.call(command)

def catandrename(filelist):
	files = []
	names = []
	with open(filelist, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			files.append(line[0])
			names.append(line[1])

	idxs = range(len(files))
	print files
	print names
	print idxs

	for idx in idxs[::2]:
		forfile1 = files[idx] + '_1.fastq.gz'
		revfile1 = files[idx] + '_2.fastq.gz'
		forfile2 = files[idx + 1] + '_1.fastq.gz'
		revfile2 = files[idx + 1] + '_2.fastq.gz'

		forfilename = names[idx] + '_1.fastq.gz'
		revfilename = names[idx] + '_2.fastq.gz'

		#Cat files

		print 'Concatentating {0} and {1} to become {2}...'.format(forfile1, forfile2, forfilename)
		command = ['cat', forfile1, forfile2]
		with open(forfilename, 'w') as outfh:
			subprocess.call(command, stdout = outfh)

		print 'Concatentating {0} and {1} to become {2}...'.format(revfile1, revfile2, revfilename)
		command = ['cat', revfile1, revfile2]
		with open(revfilename, 'w') as outfh:
			subprocess.call(command, stdout = outfh)

def renamefiles(filelist, end):
	files = []
	names = []
	with open(filelist, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			files.append(line[0])
			names.append(line[1])

	wd = '/beevol/home/taliaferro/data/TXpsi/Farris/RawReads'
	for idx, filename in enumerate(files):
		if end == 'single':
			oldfilename = os.path.join(wd, filename + '.sra.fastq.gz')
			newfilename = os.path.join(wd, names[idx] + '.fastq.gz')

			os.rename(oldfilename, newfilename)


		elif end == 'paired':
			oldfilename_for = os.path.join(wd, filename + '_1.fastq.gz')
			oldfilename_rev = os.path.join(wd, filename + '_2.fastq.gz')

			newfilename_for = os.path.join(wd, names[idx] + '_1.fastq.gz')
			newfilename_rev = os.path.join(wd, names[idx] + '_2.fastq.gz')

			os.rename(oldfilename_for, newfilename_for)
			os.rename(oldfilename_rev, newfilename_rev)


dllist(sys.argv[1], sys.argv[2])
renamefiles(sys.argv[1], sys.argv[2])
