#Usage: python Gquadruplex_cgcc.py <fasta> <windowsize> <output> <class> <region>

import re
import sys
from Bio import SeqIO
import numpy as np

def getscore(seq):
	#Given a sequence, calculate its cG/cC score.
	#First, get the longest consecutive substring of Gs and of Cs
	#http://stackoverflow.com/questions/18776238/count-the-number-of-max-consecutive-as-from-a-string-python-3

	seq = seq.upper()
	if seq.count('G') == 0:
		maxG = 0
	else:
		maxG = max(len(s) for s in re.findall(r'G+', seq))

	if seq.count('C') == 0:
		maxC = 0
	else:
		maxC = max(len(s) for s in re.findall(r'C+', seq))

	longestrun = max(maxG, maxC)

	cGscore = 0
	cCscore = 0
	#First get the cG score
	for i in range(1, longestrun + 1):
		searchstring = 'G' * i
		matches = re.findall(r'(?=({0}))'.format(searchstring), seq)
		score = len(matches) * i
		cGscore += score

	#Now the cC score
	for i in range(1, longestrun + 1):
		searchstring = 'C' * i
		matches = re.findall(r'(?=({0}))'.format(searchstring), seq)
		score = len(matches) * i
		cCscore += score

	if cCscore == 0:
		cGcCscore = cGscore
	else:
		cGcCscore = cGscore / float(cCscore)

	return cGcCscore

def iteratefasta(fasta, windowsize):
	seqcount = 0
	windowsize = int(windowsize)
	scores = []
	medianscores = {} #median score per sequence
	maxscores = {} #max score per sequence
	winners = 0
	nonwinners  = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqcount +=1
		if seqcount % 1000 == 0:
			print 'Seq {0}...'.format(seqcount)
		s = []
		seq = str(record.seq)
		if len(seq) < windowsize:
			continue
		for i in range(len(seq) - windowsize + 1):
			region = seq[i:i+windowsize]
			binnumber = round((i + 1) / float((len(seq) - windowsize + 1)), 2)
			score = getscore(region)
			scores.append([score, binnumber])
			s.append(score)

		medianscore = np.median(s)
		maxscore = max(s)
		medianscores[str(record.id).split('.')[0]] = medianscore
		maxscores[str(record.id).split('.')[0]] = maxscore
	
	return scores, medianscores, maxscores

scores, medianscores, maxscores = iteratefasta(sys.argv[1], sys.argv[2])
with open(sys.argv[3], 'a') as outfile:
	'''
	for score in scores:
		cGcCscore, binnumber = score[0], score[1]
		outfile.write(str(cGcCscore) + '\t' + str(binnumber) + '\t' + sys.argv[4] + '\t'  + sys.argv[5]  + '\n')
	'''
	outfile.write(('\t').join(['Gene', 'medianscore', 'maxscore']) + '\n')
	for record in medianscores:
		outfile.write(('\t').join([record, str(medianscores[record]), str(maxscores[record])]) + '\n')