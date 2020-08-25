import pysam
import os
import gffutils
import sys
from operator import itemgetter
from itertools import *
import numpy as np
from random import choice
import argparse

#python3

#Given a gff of regions you care about (longest 3UTR, etc.), get conservation metrics on a gene by gene basis


def getexoniccoords(gff):
	exoniccoords = {} #{gene : {chr : [[exon1start, exon1stop], [exon2start, exon2stop]]}}
	#Get the exonic coords for every gene in this gff
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)

	for gene in db.features_of_type('gene'):
		geneid = str(gene.id).replace('gene:', '').split('.')[0]
		chrm = str(gene.chrom)
		exoniccoords[geneid] = {}
		exoniccoords[geneid][chrm] = []
		for exon in db.children(gene, featuretype = 'exon', order_by = 'start'):
			exoniccoords[geneid][chrm].append([int(exon.start), int(exon.end)])

	print('Retrieved coordinates for {0} genes.'.format(len(exoniccoords)))

	return exoniccoords


def getmedianscore(coords, tbx):
	#Given a set of exonic coordinates (could be one or more exons), get the median score of all nucleotides
	#coords is the value of exoniccoords for one gene key
	#e.g. {chr2 : [[123123, 131232], [134343, 145223]]}

	scores = []
	chrm = list(coords.keys())[0]
	for exon in coords[chrm]:
		for row in tbx.fetch(chrm, exon[0], exon[1], parser = pysam.asBed()):
			score = float(row.score)
			scores.append(score)

	if not scores: #if this list is empty
		return 'NA'
	elif scores:
		medianscore = np.mean(scores)

		return medianscore

def slidingwindowmedian(coords, tbx):
	#Given a set of exonic coordinates (could be one or more exons), slide a window across the joined exons,
	#taking the median score of that window and recording it.
	#coords is the value of exoniccoords for one gene key
	#e.g. {chr2 : [[123123, 131232], [134343, 145223]]}

	medianwindowscores = []
	chrm = list(coords.keys())[0]
	joinedexoncoords = [] #every nt that is exonic, joined together
	windowsize = 100
	slidesize = 20
	for exon in coords[chrm]:
		joinedexoncoords += list(range(exon[0], exon[1] + 1))
	
	if len(joinedexoncoords) < windowsize:
		scores = []
		for coord in joinedexoncoords:
			for row in tbx.fetch(chrm, coord, coord + 1, parser = pysam.asBed()):
				scores.append(float(row.score))
		medianscore = np.mean(scores)

		return [medianscore]

	elif len(joinedexoncoords) >= windowsize:
		currentind = 0
		for ind, coord in enumerate(joinedexoncoords):
			while currentind + windowsize <= len(joinedexoncoords):
				currentwindowcoords = joinedexoncoords[currentind : currentind + windowsize]
				#print('Current window: {0}, {1}'.format(currentwindowcoords[0], currentwindowcoords[-1]))
				#OK now break this window up into chunks of consecutive integers
				#https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
				consecutivechunks = []
				for k, g in groupby(enumerate(currentwindowcoords), lambda ix : ix[0] - ix[1]):
					consecutivechunk = list(map(itemgetter(1), g))
					consecutivechunks.append([consecutivechunk[0], consecutivechunk[-1]])

				currentwindowscores = []
				for chunk in consecutivechunks:
					#print('Current chunk: {0}, {1}'.format(chunk[0], chunk[-1]))
					for row in tbx.fetch(chrm, chunk[0], chunk[-1], parser = pysam.asBed()):
						currentwindowscores.append(float(row.score))

				if len(currentwindowscores) == 0: #if there were no coords in this window that had a score
					medianwindowscores.append('NA')
				elif currentwindowscores:
					medianwindowscore = np.mean(currentwindowscores)
					medianwindowscores.append(medianwindowscore)
				currentind += slidesize

		return medianwindowscores


def iterategenes(exoniccoords, phastconsbed, genemedianfile, genewindowfile):
	#Iterate the exonic coords of each gene, getting both the median score across the whole gene 
	#and the median windowscores (there will be one for each window) for each gene

	tbx = pysam.TabixFile(phastconsbed) #phastconsbed has to have a tabix index in the same directory

	with open(genemedianfile, 'w') as medianfh, open(genewindowfile, 'w') as windowfh:
		medianfh.write(('\t').join(['ensembl_gene_id', 'medianscore']) + '\n')
		windowfh.write(('\t').join(['ensembl_gene_id', 'medianwindowscore']) + '\n')

		genecounter = 0
		for gene in exoniccoords:
			genecounter +=1
			if genecounter % 100 == 0:
				print('Analyzing gene {0} of {1}...'.format(genecounter, len(exoniccoords)))
			coords = exoniccoords[gene]
			mediangenescore = getmedianscore(coords, tbx)
			medianwindowscores = slidingwindowmedian(coords, tbx)
			medianfh.write(('\t').join([gene, str(mediangenescore)]) + '\n')
			for score in medianwindowscores:
				windowfh.write(('\t').join([gene, str(score)]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Longest feature gff.')
	parser.add_argument('--phastconsbed', type = str, help = 'Phastcons values (or phyloP values) in bed format. Can be compressed. Needs a tabix index of the file in the same directory.')
	parser.add_argument('--meanoutfile', type = str, help = 'Output file for mean scores across the entire UTR.')
	parser.add_argument('--windowoutfile', type = str, help = 'Output file for mean window scores across entire UTR.')
	args = parser.parse_args()

	exoniccoords = getexoniccoords(args.gff)
	iterategenes(exoniccoords, args.phastconsbed, args.meanoutfile, args.windowoutfile)


