#Given a gff of regions you are interested in, a fasta file of the genome sequence, and a .bed file
#of phastcons values (sorted.phastcons.mm10.bed.gz), look for kmers within the exonic regions and 
#get phastcons values from kmerstart - x to kmerend + x

from Bio import SeqIO
import sys
import pysam
import gffutils
import os
import gzip
import numpy as np

def indexgenome(genomefasta):
	sys.stderr.write('Indexing genome sequence...\n')
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	sys.stderr.write('{0} chromosomes indexed.\n'.format(len(seq_dict)))

	return seq_dict

def getSequence(seq_dict, region):
    #region = chrm;start;stop;strand
	chrm = region.split(';')[0]
	start = int(region.split(';')[1])
	stop = int(region.split(';')[2])
	strand = region.split(';')[3]

	if strand == '+':
		seq = seq_dict[chrm].seq[start:stop].upper().transcribe()
	elif strand == '-':
		seq = seq_dict[chrm].seq[start:stop].upper().reverse_complement().transcribe()

	seq = str(seq)

	return seq

def getkmerpos(featureseq, featurecoords, kmers):
	kmerindexes = []
	ntsinwindow = [] #single nt positions (in genome space) covered by a kmerofinterest
	#coords in ntinwindow are 0-based for easy interaction with a 0-based bed
	k = len(kmers[0]) #all kmers have to be the same length
	for i in range(len(featureseq) - k + 1):
		kmer = featureseq[i:i+k]
		if kmer in kmers:
			kmerindexes.append(i) # i is start of kmer; kmerofinterest goes from i to i + len(kmerofinterest)

	for ind in kmerindexes:
		for j in list(range(k)): #can add window around kmer here
			ntsinwindow.append(featurecoords[ind + j])

	return list(set(ntsinwindow))

def getscores(phastconstabix, chrom, ntsinwindow):
	#Get score for each nt in ntsinwindow
	phastconsscores = []
	for nt in ntsinwindow:
		for bed in phastconstabix.fetch(chrom, nt, nt + 1, parser = pysam.asTuple()):
			score = bed[4]
			phastconsscores.append(float(score))

	return phastconsscores

def iterategff(gff, seq_dict, kmers, windowsize, featureofinterest, phastconsbed, outfile):
	kmers = kmers.split(',')

	#Feature of interest is 2nd field of gff
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	phastconstabix = pysam.Tabixfile(phastconsbed)

	features = db.features_of_type(featureofinterest)

	with open(outfile, 'w') as outfh:
		outfh.write('Gene' + '\t' + 'medianmotifscore' + '\n')
		featurecounter = 0
		for feature in features:
			featurecounter +=1
			if featurecounter % 500 == 0:
				print('Feature {0}...'.format(featurecounter))
			featureseq = ''
			featurecoords = [] #0-based exon coords for interaction with phastcons bed
			if feature.strand == '+':
				for exon in db.children(feature, featuretype = 'exon', order_by = 'start'):
					#gff is 1-based, seq_dict is 0-based
					region = exon.chrom + ';' + str(exon.start - 1) + ';' + str(exon.end) + ';' + exon.strand
					exonseq = getSequence(seq_dict, region)
					featureseq += exonseq
					featurecoords += list(range(exon.start - 1, exon.end))

			elif feature.strand == '-':
				for exon in db.children(feature, featuretype = 'exon', order_by = 'start', reverse = True):
					#gff is 1-based
					region = exon.chrom + ';' + str(exon.start - 1) + ';' + str(exon.end) + ';' + exon.strand
					exonseq = getSequence(seq_dict, region)
					featureseq += exonseq
					featurecoords += reversed(list(range(exon.start - 1, exon.end)))

			ntsinwindow = getkmerpos(featureseq, featurecoords, kmers)
			if ntsinwindow:
				phastconsscores = getscores(phastconstabix, str(feature.chrom), ntsinwindow)
				if phastconsscores:
					score = np.median(phastconsscores)
					outfh.write(str(feature.id).split(':')[1] + '\t' + str(score) + '\n')
				elif not phastconsscores:
					outfh.write(str(feature.id).split(':')[1] + '\t' + 'NA' + '\n')
			elif not ntsinwindow:
				outfh.write(str(feature.id).split(':')[1] + '\t' + 'NA' + '\n')
			

		

seq_dict = indexgenome(sys.argv[1])
iterategff(sys.argv[2], seq_dict, 'UGUGU,GUGUG', None, 'UTR3', sys.argv[3], sys.argv[4])
#iterategff(sys.argv[2], seq_dict, 'GGGUU,UUGGG,UUUGG,GGUUU,AUUGG,GGUUA,AGGUU,UUGGA', None, 'UTR3', sys.argv[3], sys.argv[4])

