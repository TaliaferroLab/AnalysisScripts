#Window size = protein length / 40
#Step size = protein length / 100

#Usage: python leucinedensity.py CDSsequence.fasta

from Bio import SeqIO
import sys
from numpy import mean, std
import random


#Count the number of leucines int he protein.
#If each window is 1/20th of the protein, how many leucines do you expect to see in each window, assuming they are equally distributed.
def getexpectedleucinecount(AAseq, aa):
	leucinecount = AAseq.count(aa)
	proteinlength = len(AAseq)

	#Given the number of leucines in the protein, if they were equally dispersed, how many do you expect
	#to see in a given window
	expectedleucineperwindow = leucinecount / float(20)
	return leucinecount, expectedleucineperwindow

#Walk through each protein seq, counting the number of leucines in each window and comparing it to the expected number
def getrelativeleudensities(AAseq, aa):
	totalleucinecount, expectedleucineperwindow = getexpectedleucinecount(AAseq, aa)
	numberofbins = 100
	reldens = {} # {windownumber (1 to 100) : relative density}
	proteinlength = len(AAseq)
	windowsize = float(proteinlength / 20)

	#This is a little wonky.  May need to look at it later.
	stepsize = float((proteinlength - windowsize) / float(numberofbins)) * (1 + ((100/numberofbins) * 0.01))
	for i in range(numberofbins):
		windownumber = i + 1
		windowstart = int(round(i * stepsize))
		windowstop = int(round((i * stepsize) + windowsize))
		actualleucinecount = AAseq[windowstart:windowstop].count(aa)
		relativedensity = actualleucinecount / float(expectedleucineperwindow)
		reldens[windownumber] = relativedensity

	return reldens

def iteratefasta(fasta, aa):
	relativedensities = {} # {windownumber : [relative densities]}
	alldensities = []
	allcvs = [] #all cvs of subsamples for non-L aa
	for i in range(100):
		relativedensities[i + 1] = []
	for record in SeqIO.parse(fasta, 'fasta'):
		AAseq = str(record.seq.translate())[:-1] #remove '*' as stop codon
		#if there are no leucines in this sequence, skip it
		if AAseq.count(aa) == 0:
			continue
		reldens = getrelativeleudensities(AAseq, aa)
		for window in reldens:
			relativedensities[window].append(reldens[window])

	for window in relativedensities:
		#print window, mean(relativedensities[window]) #If you just want the means for a metagene
		for dens in relativedensities[window]:
			alldensities.append(dens)

	subsamplesize = len(alldensities) / 10
	for i in range(1000):
		subsampleddensities = random.sample(alldensities, subsamplesize)
		meandens = mean(subsampleddensities)
		stddens = std(subsampleddensities)
		cv = stddens/ meandens
		#if aa == 'X':
		print cv, aa, sys.argv[2]
		#else:
			#allcvs.append(cv)

	#If this isn't leucine, just print the mean of the subsampled cvs, not all of them
	#if aa != 'X':
		#print mean(allcvs), aa, sys.argv[2]


	#print 'mean = {0}, stdev = {1}, CV = {2}.'.format(mean(alldensities), std(alldensities), std(alldensities) / mean(alldensities))

	

for aa in ['G', 'P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']:
	sys.stderr.write(aa + '\n')
	iteratefasta(sys.argv[1], aa)
