import os
import subprocess
from collections import OrderedDict
import sys

def downloadfiles(readfiles):
	#Readfiles is a tab-delimited file that has the sample name and then a link to download the reads for that sample
	rfs = OrderedDict()
	with open(readfiles, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			samplename = line[0]
			fpath = line[1]
			rfs[samplename] = fpath

	samplecounter = 0
	for sample in rfs:
		samplecounter +=1
		fpath = rfs[sample]
		oldname = os.path.join(os.path.abspath('.'), fpath.split('/')[-1])
		newname = os.path.join(os.path.abspath('.'), sample + 'fastq.gz')
		print fpath, oldname
		
		print 'Downloading {0}, sample {1} of {2}...'.format(sample, samplecounter, len(rfs))
		
		#Download file
		dlcommand = ['wget', fpath]
		subprocess.call(dlcommand)

		#Rename file
		os.rename(oldname, newname)
		
	

downloadfiles(sys.argv[1])