#python3

import pysam
import os
import sys
from collections import defaultdict
import random

#Given a sam file of reads mapped to oligos and the fastq that contains the UMIs (here the first 8 nt of the reverse read)
#see how many unique UMIs there are per oligo.  Also do this after subsampling mapped reads to look at saturation.

#For a sam, get the readIDs that map to each oligoID
def oligo2readIDs(sam):
	oligo2reads = defaultdict(list) #{oligoID : [list of read IDs]}
	with pysam.AlignmentFile(sam, 'r') as infh:
		readcounter = 0
		for read in infh.fetch(until_eof = True):
			#Only consider forward reads, if single end take them all
			if read.is_paired and (read.is_read2 or read.mate_is_unmapped):
				continue
			readcounter +=1
			if readcounter % 1000000 == 0:
				print('Read {0}...'.format(readcounter))
			oligoid = read.reference_name
			readid = read.query_name
			oligo2reads[oligoid].append(readid)

	return oligo2reads

def getallreadIDs(sam):
	allreadIDs = []
	#Get all read IDs that are valid (mapped)
	with pysam.AlignmentFile(sam, 'r') as infh:
		readcounter = 0
		for read in infh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print('Read {0}...'.format(readcounter))
			#Only consider forward reads, if single end take them all
			if read.is_paired and (read.is_read2 or read.mate_is_unmapped):
				continue
			allreadIDs.append(read.query_name)

	return allreadIDs

def oligo2readIDs_sample(sam, allreadIDs, samplesize):
	samplesize = float(samplesize)
	oligo2reads = defaultdict(list) #{oligoID : [list of read IDs]}
	numbertosample = round(len(allreadIDs) * samplesize)
	print('{0} valid reads in this sam. Sampling {1}.'.format(len(allreadIDs), numbertosample))
	sampledreads = random.sample(allreadIDs, numbertosample)
	#Convert to set for faster membership tests
	sampledreads = set(sampledreads)

	#Go through the sam, only considering reads if they are in sampledreads
	with pysam.AlignmentFile(sam, 'r') as infh:
		readcounter = 0
		for read in infh.fetch(until_eof = True):
			#Only consider forward reads, if single end take them all
			if read.is_paired and (read.is_read2 or read.mate_is_unmapped):
				continue
			readcounter +=1
			if readcounter % 1000000 == 0:
				print('Read {0}...'.format(readcounter))
			readid = read.query_name
			if readid in sampledreads:
				oligoid = read.reference_name
				oligo2reads[oligoid].append(readid)

	return oligo2reads

def fastq2UMI(fastq):
	umis = {} #{read id: UMI}
	with pysam.FastxFile(fastq) as infh:
		readcounter = 0
		for entry in infh:
			readcounter +=1
			if readcounter % 1000000 == 0:
					print('Read {0}...'.format(readcounter))
			readname = entry.name
			seq = entry.sequence
			umi = seq[:8]
			umis[readname] = umi

	return umis

def umisperoligo(oligo2reads, umis, samplename):
	oligoumis = defaultdict(list) #{oligoID : [umis]}
	oligocounter = 0
	for oligo in oligo2reads:
		oligocounter +=1
		if oligocounter % 10000 == 0:
			print('Processing oligo {0} of {1}...'.format(oligocounter, len(oligo2reads)))
		readIDs = oligo2reads[oligo]
		for read in readIDs:
			readumi = umis[read]
			oligoumis[oligo].append(readumi)

	#Collapse to unique UMIs
	oligouniqueumis = {} #{oligoID : [uniqueumis]}
	for oligo in oligoumis:
		umis = oligoumis[oligo]
		uniqueumis = list(set(umis))
		oligouniqueumis[oligo] = uniqueumis

	#return oligoumis, oligouniqueumis #use this to return values for use with runsubsamples()

	
	with open(samplename + '.umis.txt', 'w') as outfh: #use this to write if you are not subsampling
		for oligo in oligo2reads:
			outfh.write(('\t').join([oligo, str(len(oligoumis[oligo])), str(len(oligouniqueumis[oligo]))]) + '\n')
	

def runsubsamples(sam, fastq, samplename):
	outfh = open(samplename + 'subsamples.txt', 'w')
	outfh.write(('\t').join(['oligo', 'totalreads', 'uniqueumis', 'samplesize']) + '\n')
	outfh.close()
	umis = fastq2UMI(fastq)
	allreadids = getallreadIDs(sam)
	for samplesize in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
		print(samplename, samplesize)
		oligo2reads = oligo2readIDs_sample(sam, allreadids, samplesize)
		oligoumis, oligouniqueumis = umisperoligo(oligo2reads, umis, samplename)
		with open(samplename + 'subsamples.txt', 'a') as outfh:
			for oligo in oligo2reads:
				outfh.write(('\t').join([oligo, str(len(oligoumis[oligo])), str(len(oligouniqueumis[oligo])), str(samplesize)]) + '\n')



#sam, fastq, outfilename
#runsubsamples(sys.argv[1], sys.argv[2], sys.argv[3])

samples = ['TWIST_HUMAN_POOL_S14_L003', 'CAF_FF_INTEGRANT_S13_L003']
samplenames = ['hg38OligoPool', 'hg38Integrants']

readdir = '/beevol/home/taliaferro/data/cisElementScreen/human/IntegrationTest/RawReads/'
samdir = '/beevol/home/taliaferro/data/cisElementScreen/human/IntegrationTest/Mapping/'

'''
#Subsampling reads
for idx, sample in enumerate(samples):
	revreads = os.path.join(readdir, sample + '_R2_001.fastq.gz')
	samplename = samplenames[idx]
	samfile = os.path.join(samdir, samplename + '.sam')
	print('Analyzing {0}, sample {1} of {2}...'.format(samplename, idx + 1, len(samplenames)))

	runsubsamples(samfile, revreads, samplename)
'''

#All reads
for idx, sample in enumerate(samples):
	revreads = os.path.join(readdir, sample + '_L004_R2_001.fastq.gz')
	samplename = samplenames[idx]
	samfile = os.path.join(samdir, samplename + '.sam')

	print('Analyzing {0}, sample {1} of {2}...'.format(samplename, idx + 1, len(samplenames)))
	print('Getting read IDs...')
	oligo2reads = oligo2readIDs(samfile)
	print('Getting UMIs...')
	umis = fastq2UMI(revreads)
	print('Connecting oligo counts and UMIs...')
	umisperoligo(oligo2reads, umis, samplename)

