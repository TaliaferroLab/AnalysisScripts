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

#Given a rnafoldGquad outfile, pull the "region" that corresponds to the G4 region (+/- x bp) and see 
#if it is more/less conserved than the rest of the UTR.

def getGquadboundaries(poslist, windowpad):
	#Rnafold gquad out files give the 1-based position of Gs that participate in quadruplex folds.
	#We have to turn this into genome coordinates.
	#So given a list of positions (e.g. [43,44,47,48,50,51,53,54,162,163,165,166,170,171,173,174]),
	#we have to split this into the two quadruplex "strcutres".
	#We'll say that if there is a break of at least 10 nt between a quadruplexed G and the previous 
	#quadruplexed G, it's the beginning of a new structure.

	structureboundaries = [] #nested list where each sublist is [start, stop] for each structure
	for idx, pos in enumerate(poslist):
		if idx == 0: #first pos
			currentstart = pos
		elif idx == len(poslist) - 1: #last pos
			currentstop = pos
			structureboundaries.append([currentstart - windowpad, currentstop + windowpad])
		elif pos - 10 > poslist[idx - 1]: #if theres a jump of more than 10 from the previous G
			currentstop = poslist[idx - 1]
			structureboundaries.append([currentstart - windowpad, currentstop + windowpad])
			currentstart = pos

	return structureboundaries

def makedb(gff):
	#Given a gff file of the region we are interested in (usually longest cds, 3' utr, etc.),
	#make a gff db so we can easily get coords
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)

	return db

def parseGquadout(gquadoutfile, db):
	#Given a RNAfold gquad out file, parse it to get gene name, feature coords, and the gquad pos
	gquaddict = {} #{genename : [chrm, strand, [exoniccoords], [[gquadregionstart, gquadregionstop]]]}

	genecounter = 0
	with open(gquadoutfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'seqname':
				continue
			genecounter +=1
			genename = line[0].split('_')[0]
			poslist = line[4].split(',')
			if poslist[0] == 'none':
				continue
			#Turn poslist into integers
			poslist = [int(pos) for pos in poslist]
			brokenposlist = getGquadboundaries(poslist, 0)
			#Get the exonic coordinates of this gene
			exoniccoords = []
			geneid = 'gene:' + genename
			chrm = db[geneid].chrom
			strand = db[geneid].strand
			for exon in db.children(geneid, featuretype = 'exon'):
				exoniccoords += range(exon.start, exon.end + 1)

			gquaddict[genename] = [chrm, strand, sorted(exoniccoords), brokenposlist]

	return gquaddict

def getGquadregions(gquaddict):
	#OK now given a dictionary where each key is a gene (that has at least 1 guad in this longest CDS, 3'UTR, whatever),
	#and the values are [chrm, strand, [exoniccoords], [[gquadregionstart, gquadregionstop]]],
	#for each gene get the exonic coordinates that lie in a gquad and those that don't.

	coords = {} #{gene : [chrm, strand, [exonic coords that are in gquad (one list for each region)], [exoniccoords that are not in gquad (one list for each uniterrupted stretch) because these can be interrupted by gquads]]}

	for gene in gquaddict:
		gquadcoords = []
		chrm = gquaddict[gene][0]
		strand = gquaddict[gene][1]
		exoniccoords = gquaddict[gene][2]
		gquadregions = gquaddict[gene][3]
		#If this region is on the minus strand, then the gquad regions start counting from the right
		if strand == '-':
			for gquadregion in gquadregions:
				gquadstart = gquadregion[1] * -1
				gquadstop = gquadregion[0] * -1
				if gquadstop == 0:
					gquadstop = -1
				gquadcoords.append(list(range(exoniccoords[gquadstart], exoniccoords[gquadstop] + 1)))

		elif strand == '+':
			for gquadregion in gquadregions:
				gquadstart = gquadregion[0]
				gquadstop = gquadregion[1]
				gquadcoords.append(list(range(exoniccoords[gquadstart], exoniccoords[gquadstop] + 1)))

		#Get a list of all nt that are in a gquad region
		allgquadnt = []
		for gquadcoord in gquadcoords:
			allgquadnt += gquadcoord

		#OK now we have the exonic coords that are in a gquad region, all other exonic coords are in non-gquad regions
		nongquadcoords = list(set(exoniccoords) - set(allgquadnt))
		nongquadcoords = sorted(nongquadcoords)
		#Now break these nongquadcoords up into chunks where the positions are consecutive
		nongquadchunks = []
		for k, g in groupby(enumerate(nongquadcoords), lambda x: x[0]-x[1]):
			nongquadchunks.append(list(map(itemgetter(1), g)))

		coords[gene] = [chrm, strand, gquadcoords, nongquadchunks]

	return coords

def getscoreofchunk(chrm, coordlist, tbx):
	#given a list of coords (a gquad chunk or a nongquad chunk) intersect with a conservation bed and get the median score in that window
	#These coordlists but be contiguous within an exon.
	phastconsvalues = []
	#Sort just in case
	coordlist = sorted(coordlist)

	for row in tbx.fetch(chrm, coordlist[0], coordlist[-1], parser = pysam.asBed()):
		score = float(row.score)
		phastconsvalues.append(score)

	medianscore = np.median(phastconsvalues)

	return medianscore

def makerandomchunks(nongquadcoordlist, chunklength, strand):
	#Given a nested list of non-gquad coords, pick a random CONTINGUOUS chunk of length chunklength

	#Get length of each contiguous exonic non g-quad region
	regionstarts = []
	regionends = []
	allpositions = [] #all possible positions to choose from
	for region in nongquadcoordlist:
		regionstarts.append(region[0])
		regionends.append(region[-1])
		#allpositions += list(range(region[0], region[1] + 1))
		allpositions += region
	#if this is on the + strand, the start of a random chunk can be anywhere in a region that is at least chunklength away from the end of that region
	trycounter = 0
	if strand == '+':
		chosen = False
		while chosen == False: #have we found a suitable random start yet?
			trycounter +=1
			#print(trycounter)
			#Pick a random position
			randomstart = choice(allpositions)
			#Make sure it isn't too close to a boundary
			tooclose = False
			for regionend in regionends:
				if randomstart + chunklength >= regionend and randomstart < regionend:
					tooclose = True

			if tooclose == False:
				chosen = True

		randomchunk = list(range(randomstart, randomstart + chunklength + 1))

	#if this is on the - strand, the end of a random chunk can be anywhere in a region that is at least chunklength away from the start of that region
	#The reason we go through this trouble is because we are going to be biased against getting things "downstream" in the transcript
	#If we didn't do this, we would have that bias for + strand genes, but would have be biased against getting "upstream" things in - strand genes
	elif strand == '-':
		chosen = False
		while chosen == False:
			trycounter +=1
			#print(trycounter)
			#Pick a random position
			randomend = choice(allpositions)
			#Make sure it isn't too close to a boundary
			tooclose = False
			for regionstart in regionstarts:
				if randomend - chunklength <= regionstart and randomend > regionstart:
					tooclose = True

			if tooclose == False:
				chosen = True

		randomchunk = list(range(randomend - chunklength, randomend + 1))

	return randomchunk

def iterategenes(coords, phastconsbed, outfile):
	#OK now that we have a dict where each key is a gene and the values contain the exonic gquadcoords and the exonic nongquadcoords,
	#make random chunks from the nongquadcoords of the same length as the gquad chunks and get conservation scores for all of them

	#coords = {} #{gene : [chrm, strand, [exonic coords that are in gquad (one list for each region)], [exoniccoords that are not in gquad (one list for each uniterrupted stretch) because these can be interrupted by gquads]]}
	scores = {} #{gene : [[gquadscore, log2ratio, percentile] (one list for each gquad in that gene)]}
	tbx = pysam.TabixFile(phastconsbed) #phastconsbed has to have a tabix index in the same directory

	with open(outfile, 'w') as outfh:
		outfh.write(('\t').join(['gene', 'gquadscore', 'medianrandomscore', 'log2scoreratio', 'gquadscorepercentile']) + '\n')

	genecounter = 0
	for gene in coords:
		genecounter +=1
		scores[gene] = []
		chrm = coords[gene][0]
		strand = coords[gene][1]
		gquadcoords = coords[gene][2]
		nongquadcoords = coords[gene][3]
		
		#Get score of gquad regions
		for chunk in gquadcoords:
			#print(gene, 'Gquadchunk', chunk)
			gquadscore = getscoreofchunk(chrm, chunk, tbx)
			randomscores = []

			#Ok now choose 500 random non-gquad chunks that were the same length as the gquad chunk
			for i in range(500):
				randomchunk = makerandomchunks(nongquadcoords, len(chunk), strand)
				#get score
				#print(i, 'randomchunk', randomchunk)
				randomscore = getscoreofchunk(chrm, randomchunk, tbx)
				randomscores.append(randomscore)

			#Now for each gquad chunk, compare it to the random chunks, getting a log2ratio between the gquad and median nongquad and an empirical pvalue
			medianrandomscore = np.median(randomscores)
			scoreratio = np.log2(gquadscore / medianrandomscore)
			#Now for pvalue, insert the gquadscore into the list of randomscores
			randomscores.append(gquadscore)
			randomscores = sorted(randomscores)
			gquadscoreindex = randomscores.index(gquadscore)
			gquadscorepercentile = gquadscoreindex / len(randomscores)
			scores[gene].append([gquadscore, scoreratio, gquadscorepercentile])
			print('Gene {0} of {1}...'.format(genecounter, len(coords)), gene, gquadscore, medianrandomscore, scoreratio, gquadscorepercentile)
			scores[gene].append([gquadscore, medianrandomscore, scoreratio, gquadscorepercentile])
			with open(outfile, 'a') as outfh:
				outfh.write(('\t').join([gene, str(gquadscore), str(medianrandomscore), str(scoreratio), str(gquadscorepercentile)]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Longest feature gff.')
	parser.add_argument('--phastconsbed', type = str, help = 'Phastcons values in bed format. Can be compressed.')
	parser.add_argument('--gquadoutfile', type = str, help = 'Output of UTRfold_gquadruplex.py.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	print('Making db...')
	db = makedb(args.gff)
	print('Parsing gquad file...')
	gquaddict = parseGquadout(args.gquadoutfile, db)
	print('Getting coords of gquads...')
	coords = getGquadregions(gquaddict)
	print('Getting conservation values...')
	iterategenes(coords, args.phastconsbed, args.outfile)
