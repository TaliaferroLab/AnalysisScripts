import os
import sys
import gzip

#Renames salmon output directories made by Austin from nonuseful IDs to IDs that contain patient and tumor IDs.
#Uses a file that relates these things (tcga_clinical_meta.tsv.gz)
#USAGE: python renamesalmon.py tcga_clinical_meta.tsv.gz <directory contianing all salmon output directories>

def getrels(idfile):
	with gzip.open(idfile, 'rt') as infh:
		rels = {} #{id : [patientid, tumor]}
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'File_ID': #skip header
				continue
			sampleid = line[0]
			tumor = line[1]
			patientid = line[2]
			rels[sampleid] = [patientid, tumor]

	return rels

def renamedirs(salmondir, rels):
	salmondirs = [os.path.join(os.path.abspath(salmondir), d) for d in os.listdir(salmondir) if os.path.isdir(os.path.join(os.path.abspath(salmondir), d)) and d.endswith('idx')]
	for salmondir in salmondirs:
		sampleid = os.path.basename(salmondir).split('_')[0]
		patientid = rels[sampleid][0]
		tumor = rels[sampleid][1]
		idx = os.path.basename(salmondir).split('_')[1]
		newpath = os.path.join(os.path.dirname(salmondir), patientid + '_' + tumor + '_' + idx)
		#Sometimes there is more than one sample per patient, or even more than one ID per sample.  Doesn't make much sense, but it's rare
		try:
			os.rename(salmondir, newpath)
		except FileExistsError:
			continue



rels = getrels(sys.argv[1])
renamedirs(sys.argv[2], rels)
