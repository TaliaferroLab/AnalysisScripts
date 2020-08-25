#If kallisto was run using transcripts supplied from a Gencode annotation, transcript IDs may have a .# on the end.
#Remove these so that transcript IDs (ENSMUST00000000) can be matched with gene IDs (ENSMUSG0000000).
#Supply the directory containing subdirectories of all samples (usually kallisto_quants).

import os
import sys

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def removedot(infile, outfile):
	with open(infile, 'r') as i, open(outfile, 'w') as o:
		for line in i:
			line = line.strip().split('\t')
			if not line[0].startswith('ENSMUST'):
				o.write(('\t').join(line) + '\n')
			else:
				tname = line[0].split('.')[0]
				o.write(('\t').join([tname, line[1], line[2], line[3], line[4]]) + '\n')

def fixfiles(quantsdir):
	for sampledir in listdir_fullpath(quantsdir):
		if os.path.isdir(sampledir):
			for f in os.listdir(sampledir):
				if os.path.basename(f) == 'quant.sf':
					infile = os.path.join(sampledir, f)
					outfile = os.path.join(sampledir, 'test.tsv')
					removedot(infile, outfile)
					os.rename(os.path.join(sampledir, 'test.tsv'), os.path.join(sampledir, 'quant.sf'))

fixfiles(sys.argv[1])