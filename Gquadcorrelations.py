#Take a fasta. Over all sequences in the fasta, get correlation of metric (delta LR, delta TPM, or delta PSI) with
#gquad density measured by four different methods.

import re
from Bio import SeqIO
import numpy as np
from scipy.stats import spearmanr, pearsonr, linregress
import sys
import argparse
import random
from itertools import izip
import math as math

def getWGGAdens(seq):
	#Given a sequence, return its WGGA-(N0-6)-WGGA-(N0-6)-WGGA-(N0-6)-WGGA dens
	seqlen = len(seq)
	motifmatches = re.findall(r'(?=([AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA))', seq)
	dens = len(motifmatches) / float(seqlen)
	return dens

def getWGGAperGGA(seq):
	ggacount = seq.count('GGA')
	motifmatches = re.findall(r'(?=([AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA))', seq)
	if ggacount == 0:
		dens = 0
	else:
		dens = len(motifmatches) / float(ggacount)
	return dens

def getcGcC(seq):
	#Do we want cGcC over the whole seq?
	#Mean over 80 bp windows?
	#Max score over all windows?
	windowsize = 80
	cGcCscores = []
	for i in range(len(seq) - windowsize + 1):
		window = seq[i:i+windowsize]

		if window.count('G') == 0:
			maxG = 0
		else:
			maxG = max(len(s) for s in re.findall(r'G+', window))

		if window.count('C') == 0:
			maxC = 0
		else:
			maxC = max(len(s) for s in re.findall(r'C+', window))

		longestrun = max(maxG, maxC)

		cGscore = 0
		cCscore = 0
		#First get the cG score
		for i in range(1, longestrun + 1):
			searchstring = 'G' * i
			matches = re.findall(r'(?=({0}))'.format(searchstring), window)
			score = len(matches) * i
			cGscore += score

		#Now the cC score
		for i in range(1, longestrun + 1):
			searchstring = 'C' * i
			matches = re.findall(r'(?=({0}))'.format(searchstring), window)
			score = len(matches) * i
			cCscore += score

		if cCscore == 0:
			cGcCscore = cGscore
		else:
			cGcCscore = cGscore / float(cCscore)

		cGcCscores.append(cGcCscore)

	meanscore = np.mean(cGcCscores)
	maxscore = max(cGcCscores)
	#return meanscore
	return maxscore

def getGCcontent(seq):
	seqlen = len(seq)
	gplusc = seq.count('G') + seq.count('C')
	gc = gplusc / float(seqlen)
	return gc

def getGoverC(seq):
	g = seq.count('G') + 1
	c = seq.count('C') + 1
	goverc = g / float(c)
	return goverc 

def getRNAfoldgquad(gquadout):
	rnafoldgquads = {} #{genename : number of G4}
	with open(gquadout, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'seqname':
				continue
			genename = line[0].split('_')[0].split('.')[0]
			seqlen = float(line[1])
			gquadpos = line[4]
			if gquadpos == 'none':
				gquads = 0
			elif ',' in gquadpos:
				gquadpos = gquadpos.split(',')
				gquads = len(gquadpos) / 8.0
			elif ',' not in gquadpos and gquadpos != 'none':
				print 'ERROR: Are you sure you have the right field for gquadpos?'
			rnafoldgquads[genename] = gquads
	return rnafoldgquads

def getdeltaLR(tpmtable):
	deltaLRs = {} #{genename : deltaLR}
	with open(tpmtable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'ensembl_gene_id':
				continue
			genename = line[0]
			deltaLR = float(line[19]) - float(line[18])
			deltaLRs[genename] = deltaLR

	return deltaLRs

def getdeltaTE(tpmtable):
	deltaTEs = {} #{genename : deltaTE}
	with open(tpmtable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'ensembl_gene_id':
				continue
			genename = line[0]
			deltaTE = float(line[49]) - float(line[47])
			deltaTEs[genename] = deltaTE

	return deltaTEs

def getserumdeltaTE(tpmtable):
	deltaTEs = {} #{genename : deltaTE}
	with open(tpmtable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'ensembl_gene_id':
				continue
			genename = line[0]
			deltaTE = float(line[47]) - float(line[46])
			deltaTEs[genename] = deltaTE

	return deltaTEs


def getSomaRatios(tpmtable):
	somaratios = {} #{genename : somaratio}
	with open(tpmtable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'ensembl_gene_id':
				continue
			genename = line[0]
			somaratio = float(line[20])
			somaratios[genename] = somaratio

	return somaratios

def getPSIvalues(psitable):
	#CADFmr1PSITable_SE.txt
	psivalues = {} #{event : deltapsi}
	with open(psitable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'Event':
				continue
			eventname = line[0]
			deltapsi = float(line[16]) - float(line[15])
			psivalues[eventname] = deltapsi
	return psivalues

def randomizedict(d):
	keys = d.keys()
	values = d.values()
	random.shuffle(keys)
	random.shuffle(values)
	shuffleddict = dict(izip(keys, values))
	return shuffleddict

def correlatescores(fasta, scoremode, metric, metrictable, gquadout):
	#This function is for looking at WGGA, WGGA per GGA, and cGcC
	#To look at RNAfoldgquad, use correlatescores_rnafoldgquad
	#fasta is fasta file of all seqs to look at
	#mode is one of ['WGGAdens', 'WGGAperGGA', 'cGcC', 'rnafold']
	#if its 'rnafold', you must supply the gquadout file, otherwise specify it as 'None'
	#metric is one of ['deltaLR', 'somaratio', 'deltapsi']
	#metrictable is tpmtable for deltaLR and somaratio and PSI table for delta PSI

	#scores are g4 prediction scores
	#metrics are either deltaLRs or somaratios
	scores = []
	metrics = []

	if metric == 'deltaLR':
		metricdict = getdeltaLR(metrictable)
	elif metric == 'somaratio':
		metricdict = getSomaRatios(metrictable)
	elif metric == 'deltapsi':
		metricdict = getPSIvalues(metrictable)
	elif metric == 'deltaTE':
		metricdict = getdeltaTE(metrictable)
	elif metric == 'serumdeltaTE':
		metricdict = getserumdeltaTE(metrictable)
	if scoremode == 'rnafold':
		rnafoldgquads = getRNAfoldgquad(gquadout)

	seqcounter = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 1000 == 0:
			print 'Sequence {0}...'.format(seqcounter)
		genename = str(record.id).split('_')[0].split('.')[0]
		seq = str(record.seq)
		if scoremode == 'WGGAdens':
			score = getWGGAdens(seq)
		elif scoremode == 'WGGAperGGA':
			score = getWGGAperGGA(seq)
		elif scoremode == 'cGcC':
			#Sequences smaller than windowsize will not have cGcC scores
			if len(seq) <= 80:
				continue
			score = getcGcC(seq)
		elif scoremode == 'rnafold':
			#Sequences smaller than 80 nt were not folded by RNAfold
			if len(seq) <= 80:
				continue
			score = rnafoldgquads[genename]
		elif scoremode == 'gc':
			score = getGCcontent(seq)
		elif scoremode == 'goverc':
			score = getGoverC(seq)
		metric = metricdict[genename]
		scores.append(score)
		metrics.append(metric)
	rho, pvalue = spearmanr(scores, metrics)[0], spearmanr(scores,metrics)[1]
	#rvalue, pvalue = linregress(scores, metrics)[2], linregress(scores, metrics)[3]
	print 'Obtained {0} scores and {1} metrics.'.format(len(scores), len(metrics))
	print 'The correlation coefficient is {0} with a p value of {1}.'.format(rho, pvalue)

def correlatescores_randomize(fasta, scoremode, metric, metrictable, gquadout):
	scores = []
	metrics = []
	validgenes = [] # list of genenames that pass length filters (usually 80 nt, see below)

	if metric == 'deltaLR':
		metricdict = getdeltaLR(metrictable)
	elif metric == 'somaratio':
		metricdict = getSomaRatios(metrictable)
	elif metric == 'deltapsi':
		metricdict = getPSIvalues(metrictable)
	elif metric == 'deltaTE':
		metricdict = getdeltaTE(metrictable)
	elif metric == 'serumdeltaTE':
		metricdict = getserumdeltaTE(metrictable)
	if scoremode == 'rnafold':
		rnafoldgquads = getRNAfoldgquad(gquadout)

	seqcounter = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 1000 == 0:
			print 'Sequence {0}...'.format(seqcounter)
		genename = str(record.id).split('_')[0].split('.')[0]
		seq = str(record.seq)
		if scoremode == 'WGGAdens':
			score = getWGGAdens(seq)
		elif scoremode == 'WGGAperGGA':
			score = getWGGAperGGA(seq)
		elif scoremode == 'cGcC':
			#Sequences smaller than windowsize will not have cGcC scores
			if len(seq) <= 80:
				continue
			score = getcGcC(seq)
		elif scoremode == 'rnafold':
			#Sequences smaller than 80 nt were not folded by RNAfold
			if len(seq) <= 80:
				continue
			score = rnafoldgquads[genename]
		elif scoremode == 'gc':
			score = getGCcontent(seq)
		elif scoremode == 'goverc':
			score = getGoverC(seq)
		metric = metricdict[genename]
		scores.append(score)
		metrics.append(metric)
		validgenes.append(genename)
	#rvalue, pvalue = spearmanr(scores, metrics)[0], spearmanr(scores,metrics)[1]
	rvalue, pvalue = linregress(scores, metrics)[2], linregress(scores, metrics)[3]
	print 'Obtained {0} scores and {1} metrics.'.format(len(scores), len(metrics))
	print 'The correlation coefficient is {0} with a p value of {1}.'.format(rvalue, pvalue)
	print '-log p = {0}'.format(math.log10(pvalue) * -1)


	#filter metricdict to remove anything that didn't get a score
	filteredmetricdict = {} #{gene : metric}
	for gene in metricdict:
		if gene in validgenes:
			filteredmetricdict[gene] = metricdict[gene]

	print len(filteredmetricdict), len(scores)

	counter = 0
	controlmetricdicts = []
	for i in range(1000):
		if i % 1000 == 0:
			print i
		shuffleddict = randomizedict(filteredmetricdict)
		controlmetricdicts.append(shuffleddict)

	controlrvalues = []
	controlpvalues = []
	for i in range(1000):
		if i % 1000 == 0:
			print 'Control {0}...'.format(i)
		#Scores can stay in the same order
		#IDs in fasta are still in the same order, but we have now randomized their metrics
		metrics = []
		randommetricdict = controlmetricdicts[i]
		for gene in validgenes:
			metric = randommetricdict[gene]
			metrics.append(metric)
		controlrvalue, controlpvalue = linregress(scores, metrics)[2], linregress(scores, metrics)[3]
		#controlrvalue, controlpvalue = spearmanr(scores, metrics)[0], spearmanr(scores, metrics)[1]
		controlrvalues.append(controlrvalue)
		controlpvalues.append(controlpvalue)

	controlmedian = np.median(controlrvalues)
	controlstd = np.std(controlrvalues)
	zscore = (rvalue - controlmedian) / controlstd
	print 'The median of the control R values was {0}. The sd was {1}.'.format(controlmedian, controlstd)
	print 'The zscore of the relationship is {0}.'.format(zscore)





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--scoremode', type = str, choices = ['WGGAdens', 'WGGAperGGA', 'cGcC', 'rnafold', 'gc', 'goverc'])
	parser.add_argument('--metric', type = str, choices = ['deltaLR', 'somaratio', 'deltapsi', 'deltaTE', 'serumdeltaTE'])
	parser.add_argument('--metrictable', type = str)
	parser.add_argument('--gquadout', type = str)
	parser.add_argument('--fasta', type = str, help = 'Fasta of sequences, usually AlldeltaLRgenes....')
	parser.add_argument('--controls', action = 'store_true', help = 'Run with randomized controls.')
	args = parser.parse_args()

	if not args.controls:
		if args.scoremode != 'rnafold':
			gquadout = None
		elif args.scoremode == 'rnafold':
			gquadout = args.gquadout

		correlatescores(args.fasta, args.scoremode, args.metric, args.metrictable, args.gquadout)

	elif args.controls:
		if args.scoremode != 'rnafold':
			gquadout = None
		elif args.scoremode == 'rnafold':
			gquadout = args.gquadout

		correlatescores_randomize(args.fasta, args.scoremode, args.metric, args.metrictable, args.gquadout)

				 

