import pysam
import os
import re
import sys
from collections import Counter
import pandas as pd
from functools import reduce

#python3
#pysam 0.15.0

def getindelpos(cigarstring):
	insertionpositions = []
	deletionpositions = []
	#Given a cigar string, return positions of insertions and deletions
	cigarstring = cigarstring.replace('M', '_M=').replace('I', '_I=').replace('D', '_D=')
	cigarstring = cigarstring.split('=')[:-1] #last one will always be empty
	#[49_M, 2_D, 5_M, 3_I, 42_M] for example
	cigarnumbers = [int(item.split('_')[0]) for item in cigarstring] #[49, 2, 5, 3, 42]
	cigarletters = [item.split('_')[1] for item in cigarstring] #[M, D, M, I, M]
	#We don't want to count any insertions when we are looking how far into the reference a particular position is
	#To do this easily, replace any number in cigarnumbers with 0 if that index in cigarletters is I
	cigarnumbers_noins = []
	for idx, number in enumerate(cigarnumbers):
		if cigarletters[idx] == 'I':
			cigarnumbers_noins.append(0)
		elif cigarletters[idx] != 'I':
			cigarnumbers_noins.append(cigarnumbers[idx])

	for idx, letter in enumerate(cigarletters):
		#Get inesrtion positions
		if letter == 'I':
			insertionlength = cigarnumbers[idx]
			previousbases = sum(cigarnumbers_noins[:idx])
			for i in range(insertionlength):
				insertionpositions.append(previousbases + i + 1)

		#Get deletion positions
		elif letter == 'D':
			deletionlength = cigarnumbers[idx]
			previousbases = sum(cigarnumbers_noins[:idx])
			for i in range(deletionlength):
				deletionpositions.append(previousbases + i + 1)

	return insertionpositions, deletionpositions


def getmismatchpos(mdzstring):
	#First find any deletions
	mdzstring = mdzstring.replace('^', '~')
	deletednt = re.findall(r"~[\D']+", mdzstring)
	for dnt in deletednt:
		mdzstring = mdzstring.replace(dnt, '=D' + str(len(dnt) - 1) + '=')

	mutas = [] #positions where there is an A where there is not supposed to be
	mutts = []
	mutcs = []
	mutgs = []

	mdzstring = mdzstring.replace('G', '=1G=').replace('C', '=1C=').replace('A', '=1A=').replace('T', '=1T=').replace('N', '=1N=')
	mdzstring = mdzstring.split('=')
	#['2', '1T', '6', '1T', '39', '1A', '13', '1A', '16', 'D4', '19', '1A', '9']
	#Remove empties that could arise from multiple mutations in a row
	mdzstring = [item for item in mdzstring if item != '']
	mdzstringnumbers = [int(item.replace('A', '').replace('T', '').replace('C', '').replace('G', '').replace('N', '').replace('D', '')) for item in mdzstring]
	#[2, 1, 6, 1, 39, 1, 13, 1, 16, 4, 19, 1, 9]

	ainds = []
	tinds = []
	cinds = []
	ginds = []
	for idx, x in enumerate(mdzstring):
		if 'A' in x:
			ainds.append(idx)
		elif 'T' in x:
			tinds.append(idx)
		elif 'C' in x:
			cinds.append(idx)
		elif 'G' in x:
			ginds.append(idx)

	for i in ainds:
		previousbases = sum(mdzstringnumbers[:i])
		mutas.append(previousbases + 1)

	for i in tinds:
		previousbases = sum(mdzstringnumbers[:i])
		mutts.append(previousbases + 1)

	for i in cinds:
		previousbases = sum(mdzstringnumbers[:i])
		mutcs.append(previousbases + 1)

	for i in ginds:
		previousbases = sum(mdzstringnumbers[:i])
		mutgs.append(previousbases + 1)

	return mutas, mutts, mutcs, mutgs

def getrecurringmuts(listofmuts_f, listofmuts_r, overlapstart, overlapend):
	recurringmuts = []
	for m in listofmuts_f:
		if m >= overlapstart and m <= overlapend:
			if m in listofmuts_r:
				recurringmuts.append(m)
		elif m < overlapstart or m > overlapend: #if it's outside the overlap we have to assume it's recurring
			recurringmuts.append(m)

	for m in listofmuts_r:
		if m >= overlapstart and m <= overlapend:
			if m in listofmuts_f:
				recurringmuts.append(m)
		elif m < overlapstart or m > overlapend:
			recurringmuts.append(m)

	#Remove duplicates
	recurringmuts = list(set(recurringmuts))

	return recurringmuts

def iteratereads(bam):
	dels = [] #positions with deletions
	ins = [] #positions with insertions
	mismatchas = [] #positions where theres an A where there shouldn't be
	mismatchts = []
	mismatchcs = []
	mismatchgs = []

	delsr = [] #total number of deletions in read
	insr = [] 
	mismatchasr = [] 
	mismatchtsr = []
	mismatchcsr = []
	mismatchgsr = []

	mismatchesperread = {} #{number of mismatches : [list of mismatches]}

	references = [] #Where are these reads mapping?

	readcounter = 0
	perfectmdz = 0
	perfectcigarstring = 0
	totallyperfectread = 0
	with pysam.AlignmentFile(bam, 'r') as infh:
		for read in infh.fetch(until_eof = True):
			if read.is_read1:
				if read.mate_is_unmapped:
					continue
				readcounter +=1
				if readcounter % 1000000 == 0:
					print('Read {0}...'.format(readcounter))
				read1start = read.reference_start + 1
				read1end = read.reference_end + 1
				cigarstring = read.cigarstring
				mdzstring = read.get_tag('MD')
				insertionpositions_r1, deletionpositions_r1 = getindelpos(cigarstring)
				mutas_r1, mutts_r1, mutcs_r1, mutgs_r1 = getmismatchpos(mdzstring)
				#These positions are relative to the READ, not the reference
				#We want things relative to the reference, so these positions must be adjusted by read.reference_start - 1
				#For most reads, this wont change anything because read.reference_start == 1
				insertionpositions_r1 = [p + (read.reference_start - 1) for p in insertionpositions_r1]
				deletionpositions_r1 = [p + (read.reference_start - 1) for p in deletionpositions_r1]
				mutas_r1 = [p + (read.reference_start - 1) for p in mutas_r1]
				mutts_r1 = [p + (read.reference_start - 1) for p in mutts_r1]
				mutcs_r1 = [p + (read.reference_start - 1) for p in mutcs_r1]
				mutgs_r1 = [p + (read.reference_start - 1) for p in mutgs_r1]

				#151bp reads, - 28 that were taken off the forward read (20nt adapter and 8 nt UMI)
				#This will have to be changed depending on the length of the insert and the lengths of the reads
				#Right now the easiest way is to look at the sam and see what the most common mdzstrings and cigarstrings are 
				#and assume that those represent perfect reads.
				if mdzstring == '104':
					perfectmdz +=1
				if cigarstring == '104M':
					perfectcigarstring +=1
				if mdzstring == '104' and cigarstring == '104M':
					totallyperfectread += 1
				references.append(read.reference_name)
			elif read.is_read2:
				if read.mate_is_unmapped:
					continue
				read2start = read.reference_start + 1
				read2end = read.reference_end + 1
				cigarstring = read.cigarstring
				mdzstring = read.get_tag('MD')
				insertionpositions_r2, deletionpositions_r2 = getindelpos(cigarstring)
				mutas_r2, mutts_r2, mutcs_r2, mutgs_r2 = getmismatchpos(mdzstring)
				#print(mutas_r2, mutts_r2, mutcs_r2, mutgs_r2)

				#All of the reverse reads mutaitons and indel positions are relative to the read, not to the reference
				#They must be adjusted ("flipped") to be relative to the reference
				#pos1 on the read = read.reference_end
				#pos2 in read = read.reference_end - 1
				#posn in read = read.reference_end - (n -1)
				insertionpositions_r2 = [read.reference_end - (p - 1) for p in insertionpositions_r2]
				deletionpositions_r2 = [read.reference_end - (p - 1) for p in deletionpositions_r2]
				mutas_r2 = [read.reference_end - (p - 1) for p in mutas_r2]
				mutts_r2 = [read.reference_end - (p - 1) for p in mutts_r2]
				mutcs_r2 = [read.reference_end - (p - 1) for p in mutcs_r2]
				mutgs_r2 = [read.reference_end - (p - 1) for p in mutgs_r2]

				#OK so if a mismatch/indel is seen as WT in the other read, don't record it as a mismatch
				#This gets tricky though because for some files the reads don't overlap all the way.
				#So for some positions you only have one observation.
				overlapstart = max(read1start, read2start)
				overlapend = min(read1end, read2end)

				#Now check if mutations that are in the overlap area are seen in both reads
				readins = getrecurringmuts(insertionpositions_r1, insertionpositions_r2, overlapstart, overlapend)
				readdels = getrecurringmuts(deletionpositions_r1, deletionpositions_r2, overlapstart, overlapend)
				readamuts = getrecurringmuts(mutas_r1, mutas_r2, overlapstart, overlapend)
				readtmuts = getrecurringmuts(mutts_r1, mutts_r2, overlapstart, overlapend)
				readcmuts = getrecurringmuts(mutcs_r1, mutcs_r2, overlapstart, overlapend)
				readgmuts = getrecurringmuts(mutgs_r1, mutgs_r2, overlapstart, overlapend)

				#Add these positions to our lists
				ins += readins
				dels += readdels
				mismatchas += readamuts
				mismatchts += readtmuts
				mismatchcs += readcmuts
				mismatchgs += readgmuts

				insr += [len(readins)]
				delsr += [len(readdels)]
				mismatchasr += [len(readamuts)]
				mismatchtsr += [len(readtmuts)]
				mismatchcsr += [len(readcmuts)]
				mismatchgsr += [len(readgmuts)]

				#total number of errors in the read
				totalerrors = []
				for error, errortype in zip([readins, readdels, readamuts, readtmuts, readcmuts, readgmuts], ['ins', 'del', 'a', 't', 'c', 'g']):
					errorlist = [errortype] * len(error) #['ins', 'ins', 'ins'] for a read with three insertions
					totalerrors += errorlist
				
				if len(totalerrors) not in mismatchesperread:
					mismatchesperread[len(totalerrors)] = []
				if len(totalerrors) == 0:
					mismatchesperread[len(totalerrors)] += ['perfect']
				else:
					mismatchesperread[len(totalerrors)] += totalerrors




	print(readcounter, perfectcigarstring, perfectmdz, totallyperfectread)
	return readcounter, references, ins, dels, mismatchas, mismatchts, mismatchcs, mismatchgs, insr, delsr, mismatchasr, mismatchgsr, mismatchcsr, mismatchtsr, mismatchesperread

def tabulateresults(readcounter, references, ins, dels, mismatchas, mismatchts, mismatchcs, mismatchgs, insr, delsr, mismatchasr, mismatchgsr, mismatchcsr, mismatchtsr, mismatchesperread):
	insdict = dict(Counter(ins))
	delsdict = dict(Counter(dels))
	adict = dict(Counter(mismatchas))
	tdict = dict(Counter(mismatchts))
	cdict = dict(Counter(mismatchcs))
	gdict = dict(Counter(mismatchgs))
	insrdict = dict(Counter(insr))
	delsrdict = dict(Counter(delsr))
	ardict = dict(Counter(mismatchasr))
	trdict = dict(Counter(mismatchtsr))
	grdict = dict(Counter(mismatchgsr))
	crdict = dict(Counter(mismatchcsr))
	referencecountdict = dict(Counter(references))

	#Get frequencies
	insdict = {key : (value / float(readcounter)) for key, value in insdict.items()}
	delsdict = {key : value / float(readcounter) for key, value in delsdict.items()}
	adict = {key : value / float(readcounter) for key, value in adict.items()}
	tdict = {key : value / float(readcounter) for key, value in tdict.items()}
	gdict = {key : value / float(readcounter) for key, value in gdict.items()}
	cdict = {key : value / float(readcounter) for key, value in cdict.items()}

	#Add missing positions and remove those over 260 (how are they there??)
	for d in [insdict, delsdict, adict, tdict, gdict, cdict]:
		d = {k : v for k, v in d.items() if k <= 260}
		for p in range(260):
			if (p + 1) not in d:
				d[p + 1] = 0.0

	print(sum(delsdict.values()))
	print(sum(insdict.values()))
	print(sum(adict.values()))
	print(sum(tdict.values()))
	print(sum(gdict.values()))
	print(sum(cdict.values()))

	#Mismatches per read df
	kdfs = []
	d = {} #{number of mismatches : {mismatch : freq}}
	for k in mismatchesperread:
		d[k] = {}
		kdict = dict(Counter(mismatchesperread[k])) #{'A' : 32, 'T': 57}
		print(k, type(k), kdict)
		muts = ['a', 't', 'g', 'c', 'ins', 'del', 'perfect']
		for mut in muts:
			if k == 0:
				print('yo')
				try:
					readnumber = kdict[mut]
					freq = readnumber / readcounter
					d[k][mut] = freq
				except KeyError:
					pass
			else:
				try:
					readnumber = kdict[mut] / k #number of reads with this type of mutation
					freq = readnumber / readcounter
					d[k][mut] = freq
				except KeyError:
					d[k][mut] = 0

	with open('MismatchPerRead.txt', 'w') as outfh:
		outfh.write(('\t').join(['numberoferrors', 'errortype', 'freq']) + '\n')
		for k in d:
			for mut in d[k]:
				outfh.write(('\t').join([str(k), mut, str(d[k][mut])]) + '\n')


	#Write frequencies per position
	insdf = pd.DataFrame.from_dict(insdict, orient = 'index', columns = ['insertionfreq'])
	deldf = pd.DataFrame.from_dict(delsdict, orient = 'index', columns = ['deletionfreq'])
	adf = pd.DataFrame.from_dict(adict, orient = 'index', columns = ['afreq'])
	tdf = pd.DataFrame.from_dict(tdict, orient = 'index', columns = ['tfreq'])
	gdf = pd.DataFrame.from_dict(gdict, orient = 'index', columns = ['gfreq'])
	cdf = pd.DataFrame.from_dict(cdict, orient = 'index', columns = ['cfreq'])
	freqdf = reduce(lambda x, y: pd.merge(x, y, left_index = True, right_index = True), [insdf, deldf, adf, tdf, gdf, cdf])
	print(freqdf.head())
	freqdf.to_csv('Frequenciesperposition.txt', sep = '\t', index = True, index_label = 'position')

	#Write number of mismatches per read
	insrdf = pd.DataFrame.from_dict(insrdict, orient = 'index', columns = ['insertioncount'])
	delsrdf = pd.DataFrame.from_dict(delsrdict, orient = 'index', columns = ['deletioncount'])
	ardf = pd.DataFrame.from_dict(ardict, orient = 'index', columns = ['acount'])
	trdf = pd.DataFrame.from_dict(trdict, orient = 'index', columns = ['tcount'])
	crdf = pd.DataFrame.from_dict(crdict, orient = 'index', columns = ['ccount'])
	grdf = pd.DataFrame.from_dict(grdict, orient = 'index', columns = ['gcount'])
	countperreaddf = reduce(lambda x, y: x.join(y, how = 'outer'), [insrdf, delsrdf, ardf, trdf, crdf, grdf])
	print(countperreaddf.head())
	countperreaddf.to_csv('Countsperread.txt', sep = '\t', index = True, index_label = 'Countperread')

	#Write distribution of sequences
	refdf = pd.DataFrame.from_dict(referencecountdict, orient = 'index', columns = ['count'])
	refdf.to_csv('ReferenceCounts.txt', sep = '\t', index = True, index_label = 'count')


readcounter, references, ins, dels, mismatchas, mismatchts, mismatchcs, mismatchgs, insr, delsr, mismatchasr, mismatchgsr, mismatchcsr, mismatchtsr, mismatchesperread = iteratereads(sys.argv[1])
tabulateresults(readcounter, references, ins, dels, mismatchas, mismatchts, mismatchcs, mismatchgs, insr, delsr, mismatchasr, mismatchgsr, mismatchcsr, mismatchtsr, mismatchesperread)
