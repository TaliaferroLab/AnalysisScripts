import pysam
from Bio import SeqIO
import sys
import pandas as pd

def countoccurrences(sam, probefasta):
	counts = {}

	#Populate counts with probenames
	with open(probefasta, 'r') as fastafh:
		for record in SeqIO.parse(fastafh, 'fasta'):
			counts[str(record.id)] = 0

	#Count how many times you see a probe in this sam.
	readcounter = 0
	unmappedcounter = 0
	with pysam.AlignmentFile(sam, 'r') as samfh:
		for read in samfh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print 'Read {0}...'.format(readcounter)
			if read.is_unmapped:
				unmappedcounter +=1
				continue
			raslprobe = read.reference_name
			counts[raslprobe] +=1

	print '{0} of {1} reads were mapped. ({2}%)'.format(readcounter - unmappedcounter, readcounter, round(((readcounter - unmappedcounter) / float(readcounter)) * 100, 2))

	#Normalize counts by the number of mapped reads
	for key, value in counts.items():
		counts[key] = value / float(readcounter - unmappedcounter)

	return counts

def iterratecountoccurrences(sams, probefasta):
	#Given a list of sam files, count the occurrences of each probe in each sam
	allcounts = {} #{sam : {probeid : counts}}
	for sam in sams:
		counts = countoccurrences(sam, probefasta)
		samname = sam.split('.')[1]
		allcounts[samname] = counts
	with open('Probecounts.txt', 'w') as outfh:
		for samname in allcounts:
			for probe in allcounts[samname]:
				counts = allcounts[samname][probe]
				outfh.write(('\t').join([samname, probe, str(counts)]) + '\n')


def getindels(sam):
	#See where indels and mismatches are in reads
	readcounter = 0
	mappedreadcounter = 0
	perfectmatches = 0
	mismatches = {} #{readname : [mismatch positions]}
	insertions = {} #{readname : [insertion positions]}
	deletions = {} #{readname : [deletion positions]}

	#This will miss reads that have both an indel and a mismatch, but that is really complicated to parse, and I'm honestly not sure I want to spend the time on that now.
	
	with pysam.AlignmentFile(sam, 'r') as samfh:
		for read in samfh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print 'Read {0}...'.format(readcounter)
			if read.is_unmapped:
				continue
			mappedreadcounter += 1
			MDtag = read.get_tag('MD')
			cigartuples = read.cigartuples
			cigaroperations = [item[0] for item in cigartuples] #0 if match, 1 if ins, 2 if deletion
			if MDtag == '46' and 1 not in cigaroperations and 2 not in cigaroperations: #this is a perfect match, first 4 nt are RYRY and were trimmed
				perfectmatches +=1
			
			#Check the cigar string for insertions and deletions
			elif 1 in cigaroperations or 2 in cigaroperations:
				if 'I' in read.cigarstring:
					insertions[read.query_name] = []
				if 'D' in read.cigarstring:
					deletions[read.query_name] = []

				currentpos = 0
				for cigartuple in cigartuples:
					operation, length = cigartuple[0], cigartuple[1]
					if operation == 0: #perfect match for this many nt
						currentpos += length
					elif operation == 1:
						for i in range(length):
							insertions[read.query_name].append(currentpos + i + 1)
						#Don't advance currentpos because these insertions aren't in the reference seq
					elif operation == 2:
						for i in range(length):
							deletions[read.query_name].append(currentpos + i + 1)
						currentpos += length

			elif MDtag != '46' and 1 not in cigaroperations and 2 not in cigaroperations:
				#then there have to be mismatches that the cigarstring is not showing, but importantly, there aren't indels here
				#so we have to look in the MD tag to see where the mismatches are
				#Turn all nucleotides into spaces so we can split on them

				mismatches[read.query_name] = []
				MDtag = MDtag.replace('A', ' ').replace('T', ' ').replace('C', ' ').replace('G', ' ').replace('N', ' ')
				MDtag = MDtag.split(' ') #e.g. [37, 3]
				MDtag = ['1' if x == '' else x for x in MDtag] #if there were two mismatches in a row
				MDtag = [int(pos) for pos in MDtag]
					
				currentpos = 0
				for ind, matchstretch in enumerate(MDtag):
					if ind + 1 < len(MDtag): #don't want to do this for the last item in MDtag
						mismatchpos = currentpos + matchstretch + 1
						mismatches[read.query_name].append(mismatchpos)
						currentpos += mismatchpos

	#print insertions, deletions, mismatches
	print 'Looked through {0} reads.  {1} were perfect matches ({2}%)'.format(mappedreadcounter, perfectmatches, round((perfectmatches / float(mappedreadcounter)) * 100, 4))
	print '{0} reads had insertions.  {1} reads had deletions.  {2} reads had mismatches.'.format(len(insertions), len(deletions), len(mismatches))

	#Reorganize dictionaries into list
	#Number of deletions, insertions, and mismatches at each position
	dels = [0] * 46
	ins = [0] * 46
	mis = [0] * 46

	#Because these were mapped against the entire oligo sequence, it is possible (although unlikely) to have indels/mismatches at position > 46.
	#This can happen if there is a deletion, such that the read carries out into reference position > 46.
	#Here, we are only going to care about the first 46 positions of the oligo.

	for readname in deletions:
		for deletion in deletions[readname]:
			if deletion > 46:
				continue
			position = deletion - 1 #make 0 based
			dels[position] +=1

	for readname in insertions:
		for insertion in insertions[readname]:
			if insertion > 46:
				continue
			position = insertion - 1
			ins[position] +=1

	for readname in mismatches:
		for mismatch in mismatches[readname]:
			if mismatch > 46:
				continue
			position = mismatch - 1
			mis[position] +=1

	positions = range(1, 47)

	#Normalize to mapped read count
	dels = [d / float(mappedreadcounter) for d in dels]
	ins = [i / float(mappedreadcounter) for i in ins]
	mis = [m / float(mappedreadcounter) for m in mis]

	#Organize into df

	df = pd.DataFrame({'position' : positions, 'insertions' : ins, 'deletions' : dels, 'mismatches' : mis, 'sample' : ['12+12cycles'] * 46})
	df.to_csv('Indels_12+12cycles.txt', sep = '\t', index = False)


#sams = sys.argv[1].split(',')
#iterratecountoccurrences(sams, sys.argv[2])

getindels(sys.argv[1])


