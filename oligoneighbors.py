#Given a set of results from the cis element screen, are the neuritelog2FC results from neighboring oligos
#more similar than results from non neighbor oligos?
#python3

import pandas as pd
import sys
from collections import defaultdict
import numpy as np
from random import shuffle

def defineoligoclasses(resultsfile):
	#For each oligo, define the oligos that are related to it in the following classes:
	#1. Neighbor (large overlap)
	#2. NextToNeighbor (small overlap)
	#3. Same gene non overlap
	#4. Different gene
	log2FC = {} #{oligoid : neuritelog2FC}
	relationships = {} #{oligoid : [[list of deltalog2FC in class 1], [list of deltalog2FC in class 2], [list of deltalog2FC in class 3], [list of deltalog2FC in class 4]]}
	signs = {} #{oligoid : [[fraction of oligos in class 1 that have same sign of log2FC], [class2], [class3], [class4]]}
	allcomparisons = [] #all x vs y, do they have the same sign of log2FC? [True, True, False, False, etc.]

	oligos = [] #all oligos in resultsfile
	with open(resultsfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'oligo':
				continue
			oligos.append(line[0])
			log2FC[line[0]] = float(line[8])

	oligocounter = 0
	for oligox in oligos:
		oligocounter +=1
		if oligocounter % 100 == 0:
			print('Oligo {0}...'.format(oligocounter))
		relationships[oligox] = [[], [], [], []]
		signs[oligox] = [[], [], [], []]
		genenamex = oligox.split('.')[0]
		oligonumberx = int(oligox.split('.')[1].split('|')[0])
		log2FCx = log2FC[oligox]
		for oligoy in oligos:
			if oligox == oligoy:
				continue
			genenamey = oligoy.split('.')[0]
			oligonumbery = int(oligoy.split('.')[1].split('|')[0])
			log2FCy = log2FC[oligoy]

			deltalog2FC = abs(log2FCy - log2FCx)
			if (log2FCx > 0 and log2FCy >0) or (log2FCx < 0 and log2FCy < 0):
				samesign = True
			else:
				samesign = False
			allcomparisons.append(samesign)
			if genenamex != genenamey: #class 4
				relationships[oligox][3].append(deltalog2FC)
				signs[oligox][3].append(samesign)
			elif genenamex == genenamey:
				if abs(oligonumberx - oligonumbery) == 1: #class 1
					relationships[oligox][0].append(deltalog2FC)
					signs[oligox][0].append(samesign)
				elif abs(oligonumberx - oligonumbery) == 2: #class 2
					relationships[oligox][1].append(deltalog2FC)
					signs[oligox][1].append(samesign)
				elif abs(oligonumberx - oligonumbery) > 2: #class 3
					relationships[oligox][2].append(deltalog2FC)
					signs[oligox][2].append(samesign)

	with open('test.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'neighbormedian', 'nexttoneighbormedian', 'samegenemedian', 'differentgenemedian']) + '\n')
		for oligo in relationships:
			class1 = np.median(relationships[oligo][0])
			class2 = np.median(relationships[oligo][1])
			class3 = np.median(relationships[oligo][2])
			class4 = np.median(relationships[oligo][3])
			outfh.write(('\t').join([oligo, str(class1), str(class2), str(class3), str(class4)]) + '\n')

	
	with open('test2.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'neighbormedian', 'nexttoneighbormedian', 'samegenemedian', 'differentgenemedian']) + '\n')
		for oligo in relationships:
			if signs[oligo][0] and signs[oligo][1] and signs[oligo][2] and signs[oligo][3]: #exclude oligos that have 0 oligos in at least 1 of the 4 classes
				class1 = sum(signs[oligo][0]) / len(signs[oligo][0])
				class2 = sum(signs[oligo][1]) / len(signs[oligo][1])
				class3 = sum(signs[oligo][2]) / len(signs[oligo][2])
				class4 = sum(signs[oligo][3]) / len(signs[oligo][3])
				outfh.write(('\t').join([oligo, str(class1), str(class2), str(class3), str(class4)]) + '\n')

	allc = sum(allcomparisons) / len(allcomparisons)
	print(allc)
	

	return relationships

def getcompleteoligos(resultsfile):
	#Get the genes for which we have a recorded log2FC for every oligo (or at least those with no gaps)
	numbers = defaultdict(list) #{oligo : [numbers that we have data for]}
	completeoligos = {} #oligos for which we have no gaps, {oligo : their log2FC IN ORDER}
	log2fcs = defaultdict(dict) #{oligoname : {number : log2fc}}
	df = pd.read_csv(resultsfile, sep = '\t', header = 0)
	for index, row in df.iterrows():
		oligo = row['oligo']
		oligoname = oligo.split('.')[0]
		oligonumber = int(oligo.split('.')[1].split('|')[0])
		log2fc = row['neuritelog2FC']
		log2fcs[oligoname][oligonumber] = log2fc

	#Go through and see if there are any breaks in oligos
	for oligo in log2fcs:
		oligonumbers = sorted(log2fcs[oligo].keys())
		if len(oligonumbers) - 1 == max(oligonumbers) - min(oligonumbers):
			completeoligos[oligo] = []
			for number in oligonumbers:
				completeoligos[oligo].append(log2fcs[oligo][number])

	#print(completeoligos)
	#print(len(completeoligos))
	return completeoligos

def getneighbordistances(l):
	#given a list of log2FC in order across a gene, get the difference in log2FC between each neighbor
	differences = [t - s for s, t in zip(l, l[1:])]
	return differences

def getzs(completeoligos):
	zs = {} #{oligo : zscore}
	oligocounter = 0
	for oligo in completeoligos:
		oligocounter +=1
		if oligocounter % 1000 == 0:
			print('Oligo {0}...'.format(oligocounter))
		log2fcs = completeoligos[oligo]
		realdistances = getneighbordistances(log2fcs)
		medianreal = np.median(realdistances)
		#shuffle 500 times
		randomdistances = []
		for i in range(500):
			shuffle(log2fcs)
			randomdistance = getneighbordistances(log2fcs)
			medrandom = np.median(randomdistance)
			randomdistances.append(medrandom)
		meanrandom = np.mean(randomdistances)
		sdrandom = np.std(randomdistances)
		z = (medianreal - meanrandom) / sdrandom
		zs[oligo] = z

	with open('zs.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'z']) + '\n')
		for oligo in zs:
			outfh.write(('\t').join([oligo, str(zs[oligo])]) + '\n')


def getneuriteenrichedneighbors(resultsfile):
	#Go through and get neighbors and nexttoneighbors of sig enriched genes
	neighbors = {} #{oligoid : class} #classes are neuriteneighbor, somaneighbor, neuritenexttoneighbor, somanexttoneighbor, none.
	#They have that listed priority, meaning that if something is both a neuriteneighbor and a somaneighbor, it will be stored here as a neuriteneighbor

	#Get sig oligos
	df = pd.read_csv(resultsfile, sep = '\t', header = 0)

	neuriteoligos = df[df['sig'] == 'neurite']['oligo'].tolist()
	somaoligos = df[df['sig'] == 'soma']['oligo'].tolist()

	for index, row in df.iterrows():
		oligo = row['oligo']
		if row['sig'] == 'neurite' or row['sig'] == 'soma':
			neighbors[oligo] = row['sig']
			continue
		oligoname = oligo.split('.')[0]
		oligonumber = int(oligo.split('.')[1].split('|')[0])
		oligogene = oligo.split('|')[1]
		oligoneighbors = [oligoname + '.' + str(oligonumber - 1) + '|' + oligogene, oligoname + '.' + str(oligonumber + 1) + '|' + oligogene]
		oligonexttoneighbors = [oligoname + '.' + str(oligonumber - 2) + '|' + oligogene, oligoname + '.' + str(oligonumber + 2) + '|' + oligogene]

		#class1 = neuriteneighbor, class2 = somaneighbor, class3 = neuritenexttoneighbor, class4 = somanexttoneighbor, class5 = none
		oligoclass = [5]
		for o in oligoneighbors:
			if o in neuriteoligos:
				oligoclass.append(1)
			if o in somaoligos:
				oligoclass.append(2)
		for o in oligonexttoneighbors:
			if o in neuriteoligos:
				oligoclass.append(3)
			if o in somaoligos:
				oligoclass.append(4)

		oligoclass = min(oligoclass)
		if oligoclass == 1:
			neighbors[oligo] = 'neuriteneighbor'
		elif oligoclass == 2:
			neighbors[oligo] = 'somaneighbor'
		elif oligoclass == 3:
			neighbors[oligo] = 'neuritenexttoneighbor'
		elif oligoclass == 4:
			neighbors[oligo] = 'somanexttoneighbor'
		elif oligoclass == 5:
			neighbors[oligo] = 'none'

	with open('neighbors.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'class']) + '\n')
		for oligo in neighbors:
			outfh.write(('\t').join([oligo, neighbors[oligo]]) + '\n')

#For every oligo, pick 2 other oligos at random.  Assess whether one or both has the same sign of log2FC.
def randomoligosigns(resultsfile):
	randompairs_diffgenes = {} #{oligo : [1.0, 0.0, 0.5, etc.]}, ten values per oligo
	randompairs_samegenes = {} #{oligo : [1.0, 0.0, 0.5, etc.]}, ten values per oligo

	df = pd.read_csv(resultsfile, sep = '\t', header = 0)
	df[['oligo', 'Gene']] = df.oligo.str.split('|', expand = True)

	oligocounter = 0
	#Only look at oligo NOT in this oligo's gene
	for index, row in df.iterrows():
		oligocounter +=1
		if oligocounter % 1000 == 0:
			print('Oligo {0}...'.format(oligocounter))
		oligo = row['oligo']
		gene = row['Gene']
		randompairs_diffgenes[oligo + '|' + gene] = []
		log2FC = row['neuritelog2FC']
		newdf = df[df['Gene'] != gene]
		if log2FC > 0:
			oligosign = 'p'
		else:
			oligosign = 'n'
		for i in range(10):
			sampdf = newdf.sample(n = 2)
			samplog2fc = sampdf['neuritelog2FC'].tolist()
			randomsigns = ['p' if s > 0 else 'n' for s in samplog2fc]
			if randomsigns.count('p') == 1:
				oligovalue = 0.5
			elif randomsigns.count('p') == 2 and oligosign == 'p':
				oligovalue = 1.0
			elif randomsigns.count('n') == 2 and oligosign == 'n':
				oligovalue = 1.0
			elif randomsigns.count('p') == 2 and oligosign == 'n':
				oligovalue = 0.0
			elif randomsigns.count('n') == 2 and oligosign == 'p':
				oligovalue = 0.0
			randompairs_diffgenes[oligo + '|' + gene].append(oligovalue)

	with open('randompairs_diffgenes.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'value']) + '\n')
		for oligo in randompairs_diffgenes:
			for value in randompairs_diffgenes[oligo]:
				outfh.write(('\t').join([oligo, str(value)]) + '\n')

	oligocounter = 0
	#Only look at oligo NOT in this oligo's gene
	for index, row in df.iterrows():
		oligocounter +=1
		if oligocounter % 1000 == 0:
			print('Oligo {0}...'.format(oligocounter))
		oligo = row['oligo']
		gene = row['Gene']
		randompairs_samegenes[oligo + '|' + gene] = []
		log2FC = row['neuritelog2FC']
		newdf = df[df['Gene'] == gene]
		if log2FC > 0:
			oligosign = 'p'
		else:
			oligosign = 'n'
		for i in range(10):
			sampdf = newdf.sample(n = 2, replace = True)
			samplog2fc = sampdf['neuritelog2FC'].tolist()
			randomsigns = ['p' if s > 0 else 'n' for s in samplog2fc]
			if randomsigns.count('p') == 1:
				oligovalue = 0.5
			elif randomsigns.count('p') == 2 and oligosign == 'p':
				oligovalue = 1.0
			elif randomsigns.count('n') == 2 and oligosign == 'n':
				oligovalue = 1.0
			elif randomsigns.count('p') == 2 and oligosign == 'n':
				oligovalue = 0.0
			elif randomsigns.count('n') == 2 and oligosign == 'p':
				oligovalue = 0.0
			randompairs_samegenes[oligo + '|' + gene].append(oligovalue)

	with open('randompairs_samegenes.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'value']) + '\n')
		for oligo in randompairs_samegenes:
			for value in randompairs_samegenes[oligo]:
				outfh.write(('\t').join([oligo, str(value)]) + '\n')






#defineoligoclasses(sys.argv[1])

#completeoligos = getcompleteoligos(sys.argv[1])
#getzs(completeoligos)

#getneuriteenrichedneighbors(sys.argv[1])
randomoligosigns(sys.argv[1])