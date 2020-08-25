#With this script, we want to take an open reading frame, then try as hard as we can to make a stretch
#within that open reading frame that has extremely high or extremely low GC content, but without changing
#any amino acids.

#python 3

import sys
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import random
from collections import OrderedDict
from numpy.random import choice

def makecodoncode():
	bases = 'TCAG'
	codons = [a + b + c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))
	return codon_table

def translateseq(seq, codon_table):
	aaseq = ''
	for i in range(0, len(seq), 3):
		codon = seq[i : i + 3]
		aa = codon_table[codon]
		aaseq += aa

	return aaseq

def getGCwindows(seq, windowsize):
	gcdict = OrderedDict() #{windowstart : gccontent of window}
	for i in range(len(seq)):
		if i + windowsize < len(seq):
			window = seq[i : i + windowsize]
			windowgc = round((window.count('C') + window.count('G')) / len(window), 2)
			gcdict[i] = windowgc

	return gcdict

def getGCwindows_centered(seq, windowsize):
	#Get GC of window centered on position
	gcdict = OrderedDict()
	for i in range(len(seq)):
		windowstart = i - (windowsize / 2)
		windowend = i + (windowsize / 2)
		if windowstart < 0:
			windowstart = 0
		if windowend > len(seq) - 1:
			windowend = len(seq) - 1

		window = seq[int(windowstart) : int(windowend)]
		windowgc = round(window.count('C') + window.count('G') / len(window), 2)
		gcdict[i] = windowgc

	return gcdict

def getmaxGCcodon(aa, codon_table):
	possiblecodons = []
	for codon in codon_table:
		if codon_table[codon] == aa:
			possiblecodons.append(codon)

	gcs = {} #{gc : [codons]}
	for codon in possiblecodons:
		gc = codon.count('G') + codon.count('C')
		if gc not in gcs:
			gcs[gc] = [codon]
		else:
			gcs[gc].append(codon)

	maxgc = max(gcs.keys())
	maxgccodons = gcs[maxgc]
	maxgccodon = random.choice(maxgccodons)

	return maxgccodon

def getminGCcodon(aa, codon_table):
	possiblecodons = []
	for codon in codon_table:
		if codon_table[codon] == aa:
			possiblecodons.append(codon)

	gcs = {} #{gc : [codons]}
	for codon in possiblecodons:
		gc = codon.count('G') + codon.count('C')
		if gc not in gcs:
			gcs[gc] = [codon]
		else:
			gcs[gc].append(codon)

	mingc = min(gcs.keys())
	mingccodons = gcs[mingc]
	mingccodon = random.choice(mingccodons)

	return mingccodon

#Initiation buffer 
def makeminseq(orfseq, codon_table, initiationbuffer):
	newseq = ''
	for i in range(0, len(orfseq), 3):
		codon = orfseq[i : i + 3]
		if i <= initiationbuffer:
			newseq += codon
		else:
			aa = codon_table[codon]
			mingccodon = getminGCcodon(aa, codon_table)
			newseq += mingccodon

	return newseq

def makemaxseq(orfseq, codon_table, initiationbuffer):
	newseq = ''
	for i in range(0, len(orfseq), 3):
		codon = orfseq[i : i + 3]
		if i <= initiationbuffer:
			newseq += codon
		else:
			aa = codon_table[codon]
			maxgccodon = getmaxGCcodon(aa, codon_table)
			newseq += maxgccodon

	return newseq


#From this code, it looks like we can make a particularly GCpoor region centered at nt 400.
#I want to make a GC poor region that is 300 nt (100 codons) wide, centered here (goes from 250 to 550), but that has the rest of the ORF untouched.
#Then I want to make a control that has the same overall GC content as the version with this window, but with the GC change more spread out.
def iscodonchangeable(orfseq, codon_table, minormax):
	#Is it possible to make a new version of this codon that has one less G/C?
	#Return codon positions (0-based) of codons where this is possible
	changeablecodons = []
	for i in range(0, len(orfseq), 3):
		changeable = False
		codon = orfseq[i : i + 3]
		codongc = codon.count('G') + codon.count('C')
		aa = codon_table[codon]
		
		possiblecodons = []
		for codon in codon_table:
			if codon_table[codon] == aa:
				possiblecodons.append(codon)

		for possiblecodon in possiblecodons:
			possiblecodongc = possiblecodon.count('G') + possiblecodon.count('C')
			if minormax == 'min':
				if possiblecodongc + 1 == codongc:
					changeable = True
			elif minormax == 'max':
				if possiblecodongc - 1 == codongc:
					changeable = True

		if changeable == True:
			changeablecodons.append(i / 3)

	return changeablecodons

#Make a version that has a GC poor/rich window, as well as a control that has some overall GC content as rich/poor window version
def makeGCpoorwindowseq(orfseq, codontable, initiationbuffer, minormax):
	changeablecodons = iscodonchangeable(orfseq, codontable, minormax)

	newseq = ''
	for i in range(0, len(orfseq), 3):
		codon = orfseq[i : i + 3]
		#if i < 1174 or i >= 1324:
		#if i < 724 or i >= 874:
		#if i < 649 or i >= 949:
		if i < 249 or i >= 549: #100 codons wide
			newseq += codon
		else:
			aa = codon_table[codon]
			if minormax == 'min':
				mingccodon = getminGCcodon(aa, codon_table)
				newseq += mingccodon
			elif minormax == 'max':
				maxgccodon = getmaxGCcodon(aa, codon_table)
				newseq += maxgccodon

	
	#Make control window seq
	controlnewseqgc = 0
	newseqgc = newseq.count('C') + newseq.count('G')
	oldseqgc = orfseq.count('C') + orfseq.count('G')

	gcdiff = abs(oldseqgc - newseqgc)

	numberofcodonswecanchange = len(changeablecodons)
	probofchange = gcdiff / numberofcodonswecanchange
	choices = ['change', 'nochange']
	weights = [probofchange, 1 - probofchange]

	#Try until you get a controlnewseq with the same GC as newseq
	while controlnewseqgc != newseqgc:
		controlnewseq = ''
		for i in range(0, len(orfseq), 3):
			codon = orfseq[i : i + 3]
			codonnumber = i / 3
			if i <= initiationbuffer:
				controlnewseq += codon
			else:
				if codonnumber in changeablecodons:
					shouldchange = choice(choices, p = weights)
					if shouldchange == 'nochange':
						controlnewseq += codon
					elif shouldchange == 'change':
						aa = codon_table[codon]
			
						possiblecodons = []
						for pcodon in codon_table:
							if codon_table[pcodon] == aa:
								possiblecodons.append(pcodon)

						possiblegccodons = []
						for possiblecodon in possiblecodons:
							codongc = codon.count('G') + codon.count('C')
							possiblecodongc = possiblecodon.count('G') + possiblecodon.count('C')
							if minormax == 'min':
								if possiblecodongc + 1 == codongc:
									possiblegccodons.append(possiblecodon)
							elif minormax == 'max':
								if possiblecodongc - 1 == codongc:
									possiblegccodons.append(possiblecodon)
						codontoswitchto = random.choice(possiblegccodons)
						controlnewseq += codontoswitchto

				elif codonnumber not in changeablecodons:
					controlnewseq += codon

		controlnewseqgc = controlnewseq.count('G') + controlnewseq.count('C')

	return newseq, controlnewseq




#For plotting changes of GC content across the entire ORF
codon_table = makecodoncode()

fireflyseq = 'atggaagacgccaaaaacataaagaaaggcccggcgccattctatccgctggaagatggaaccgctggagagcaactgcataaggctatgaagagatacgccctggttcctggaacaattgcttttacagatgcacatatcgaggtggacatcacttacgctgagtacttcgaaatgtccgttcggttggcagaagctatgaaacgatatgggctgaatacaaatcacagaatcgtcgtatgcagtgaaaactctcttcaattctttatgccggtgttgggcgcgttatttatcggagttgcagttgcgcccgcgaacgacatttataatgaacgtgaattgctcaacagtatgggcatttcgcagcctaccgtggtgttcgtttccaaaaaggggttgcaaaaaattttgaacgtgcaaaaaaagctcccaatcatccaaaaaattattatcatggattctaaaacggattaccagggatttcagtcgatgtacacgttcgtcacatctcatctacctcccggttttaatgaatacgattttgtgccagagtccttcgatagggacaagacaattgcactgatcatgaactcctctggatctactggtctgcctaaaggtgtcgctctgcctcatagaactgcctgcgtgagattctcgcatgccagagatcctatttttggcaatcaaatcattccggatactgcgattttaagtgttgttccattccatcacggttttggaatgtttactacactcggatatttgatatgtggatttcgagtcgtcttaatgtatagatttgaagaagagctgtttctgaggagccttcaggattacaagattcaaagtgcgctgctggtgccaaccctattctccttcttcgccaaaagcactctgattgacaaatacgatttatctaatttacacgaaattgcttctggtggcgctcccctctctaaggaagtcggggaagcggttgccaagaggttccatctgccaggtatcaggcaaggatatgggctcactgagactacatcagctattctgattacacccgagggggatgataaaccgggcgcggtcggtaaagttgttccattttttgaagcgaaggttgtggatctggataccgggaaaacgctgggcgttaatcaaagaggcgaactgtgtgtgagaggtcctatgattatgtccggttatgtaaacaatccggaagcgaccaacgccttgattgacaaggatggatggctacattctggagacatagcttactgggacgaagacgaacacttcttcatcgttgaccgcctgaagtctctgattaagtacaaaggctatcaggtggctcccgctgaattggaatccatcttgctccaacaccccaacatcttcgacgcaggtgtcgcaggtcttcccgacgatgacgccggtgaacttcccgccgccgttgttgttttggagcacggaaagacgatgacggaaaaagagatcgtggattacgtcgccagtcaagtaacaaccgcgaaaaagttgcgcggaggagttgtgtttgtggacgaagtaccgaaaggtcttaccggaaaactcgacgcaagaaaaatcagagagatcctcataaaggccaagaagggcggaaagatcgccgtgtaa'.upper()
'''
minfireflyseq = makeminseq(fireflyseq, codon_table, 180)
maxfireflyseq = makemaxseq(fireflyseq, codon_table, 180)

wtgcdict = getGCwindows_centered(fireflyseq, 150)
mingcdict = getGCwindows_centered(minfireflyseq, 150)
maxgcdict = getGCwindows_centered(maxfireflyseq, 150)

poscount = len(list(wtgcdict.keys()))
d = {'pos' : list(wtgcdict.keys()) * 3, 'gc' : list(wtgcdict.values()) + list(mingcdict.values()) + list(maxgcdict.values()), 'samp' : ['wt'] * poscount + ['min'] * poscount + ['max'] * poscount}

df = pd.DataFrame.from_dict(d)

ax = sns.lineplot(x = 'pos', y = 'gc', hue = 'samp', data = d)
plt.show()
'''


#For plotting changes of GC content in windows
codon_table = makecodoncode()

fireflyseq = 'atggaagacgccaaaaacataaagaaaggcccggcgccattctatccgctggaagatggaaccgctggagagcaactgcataaggctatgaagagatacgccctggttcctggaacaattgcttttacagatgcacatatcgaggtggacatcacttacgctgagtacttcgaaatgtccgttcggttggcagaagctatgaaacgatatgggctgaatacaaatcacagaatcgtcgtatgcagtgaaaactctcttcaattctttatgccggtgttgggcgcgttatttatcggagttgcagttgcgcccgcgaacgacatttataatgaacgtgaattgctcaacagtatgggcatttcgcagcctaccgtggtgttcgtttccaaaaaggggttgcaaaaaattttgaacgtgcaaaaaaagctcccaatcatccaaaaaattattatcatggattctaaaacggattaccagggatttcagtcgatgtacacgttcgtcacatctcatctacctcccggttttaatgaatacgattttgtgccagagtccttcgatagggacaagacaattgcactgatcatgaactcctctggatctactggtctgcctaaaggtgtcgctctgcctcatagaactgcctgcgtgagattctcgcatgccagagatcctatttttggcaatcaaatcattccggatactgcgattttaagtgttgttccattccatcacggttttggaatgtttactacactcggatatttgatatgtggatttcgagtcgtcttaatgtatagatttgaagaagagctgtttctgaggagccttcaggattacaagattcaaagtgcgctgctggtgccaaccctattctccttcttcgccaaaagcactctgattgacaaatacgatttatctaatttacacgaaattgcttctggtggcgctcccctctctaaggaagtcggggaagcggttgccaagaggttccatctgccaggtatcaggcaaggatatgggctcactgagactacatcagctattctgattacacccgagggggatgataaaccgggcgcggtcggtaaagttgttccattttttgaagcgaaggttgtggatctggataccgggaaaacgctgggcgttaatcaaagaggcgaactgtgtgtgagaggtcctatgattatgtccggttatgtaaacaatccggaagcgaccaacgccttgattgacaaggatggatggctacattctggagacatagcttactgggacgaagacgaacacttcttcatcgttgaccgcctgaagtctctgattaagtacaaaggctatcaggtggctcccgctgaattggaatccatcttgctccaacaccccaacatcttcgacgcaggtgtcgcaggtcttcccgacgatgacgccggtgaacttcccgccgccgttgttgttttggagcacggaaagacgatgacggaaaaagagatcgtggattacgtcgccagtcaagtaacaaccgcgaaaaagttgcgcggaggagttgtgtttgtggacgaagtaccgaaaggtcttaccggaaaactcgacgcaagaaaaatcagagagatcctcataaaggccaagaagggcggaaagatcgccgtgtaa'.upper()
newseq, controlnewseq = makeGCpoorwindowseq(fireflyseq, codon_table, 180, 'max')

fireflyaa = translateseq(fireflyseq, codon_table)
newseqaa = translateseq(newseq, codon_table)
controlseqaa = translateseq(controlnewseq, codon_table)

print(fireflyaa == newseqaa == controlseqaa)

wtgcdict = getGCwindows_centered(fireflyseq, 150)
newgcdict = getGCwindows_centered(newseq, 150)
controlgcdict = getGCwindows_centered(controlnewseq, 150)

poscount = len(list(wtgcdict.keys()))
d = {'pos' : list(wtgcdict.keys()) * 3, 'gc' : list(wtgcdict.values()) + list(newgcdict.values()) + list(controlgcdict.values()), 'samp' : ['wt'] * poscount + ['new'] * poscount + ['ctrl'] * poscount}

df = pd.DataFrame.from_dict(d)
ax = sns.lineplot(x = 'pos', y = 'gc', hue = 'samp', data = d)
plt.show()


