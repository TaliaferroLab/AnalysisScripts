#Split a fastq in n equal parts.
import os
import subprocess
import sys

def unzipfile(f):
	command = ['gunzip', f]
	subprocess.call(command)

def zipfile(f):
	command = ['gzip', f]
	subprocess.call(command)

def getnumberoflines(f):
	with open(f, 'r') as infh:
		for i, l in enumerate(infh):
			pass
		return i + 1

def splitfiles(d, parts):
	parts = int(parts)
	fcounter = 0
	fs = os.listdir(d)
	#only gzipped fastq files
	fs = [f for f in fs if f.endswith('.fastq.gz')]
	for f in fs:
		fcounter +=1
		print 'Splitting file {0}, number {1} of {2}...'.format(f, fcounter, len(fs))
		print 'Unzipping...'
		unzipfile(f)
		unzippedfile = os.path.basename(f).split('.')
		unzippedfile = ('.').join(unzippedfile[:-1])
		print 'Splitting...'
		readnumber = int(getnumberoflines(unzippedfile) / 4.0)
		readsperpart = int(int(readnumber + parts - 1) / float(parts))
		linesperpart = readsperpart * 4
		print '{0} reads in this file. Splitting to give {1} reads per file.'.format(readnumber, readsperpart)

		command = ['split', '-l', str(linesperpart), unzippedfile, unzippedfile]
		subprocess.call(command)

		for suffix in ['aa', 'ab', 'ac', 'ad']:
			fname = unzippedfile + suffix
			zipfile(fname)

		zipfile(unzippedfile)

splitfiles(sys.argv[1], sys.argv[2])