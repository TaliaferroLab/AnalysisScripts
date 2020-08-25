#Make a metaplot of reads across UTRs and CDS.  For this, it's best to use only one CDS or UTR per
#gene, otherwise you won't know what bin to assign the read to.  Going to use the longest CDS/UTR
#for each gene.

import gffutils
import os
from Bio import SeqIO
import argparse
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import gzip
import sys
import pysam
import cPickle as pickle
import math
import argparse


#Given a bam, get the ribosome "positions" (15 nt offset from read start)
#Importantly, due to introns, this may not just be readstart + 15!
def getribopositions(bam):
	readcounter = 0
	consideredreadcounter = 0
	ribopositions = {} #{chrm : {strand : [positions]}}
	if 'tot' in bam:
		acceptedreadlengths = [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
	elif 'rpf' in bam:
		acceptedreadlengths = [26, 27, 28, 29, 30, 31]
	elif 'RFP' in bam:
		print 'hello'
		acceptedreadlengths = [30, 31, 32, 33, 34, 35]
	elif 'rnaseq' in bam:
		acceptedreadlengths = range(32, 51)
	else: #for test.bam
		acceptedreadlengths = [26, 27, 28, 29, 30, 31]
	with pysam.AlignmentFile(bam, 'rb') as infh:
		for read in infh.fetch(until_eof = True):
			#Only consider reads of certain lengths
			#First read of mate pair is on the wrong strand
			if read.query_length not in acceptedreadlengths or read.is_read1:
				continue

			#Only consider unique mappers
			numberofalignments = read.get_tag('NH')
			if numberofalignments > 1 or read.is_secondary:
				continue

			readcounter +=1

			chrm = read.reference_name
			if read.is_reverse:
				strand = '-'
			elif not read.is_reverse:
				strand = '+'
			
			#We are going to get an ordered list of reference positions that this read maps to.
			#This is done with get_reference_positions()
			refpositions = read.get_reference_positions(full_length = True)
			#Using full_length = True will allow soft-clipped nts to have placeholders of None in this list
			#We still want them as placeholders so that the offset will be correct
			refpositions = [pos + 1 if type(pos) == int else pos for pos in refpositions] #make 1 based
			if read.query_length == 31 or read.query_length == 32 or read.query_length == 30 or read.query_length == 33 or read.query_length == 34 or read.query_length == 35:
				if strand == '+':
					ribosomeposition = refpositions[13]
				elif strand == '-':
					refpositions.reverse()
					ribosomeposition = refpositions[13]
			else:
				if strand == '+':
					ribosomeposition = refpositions[12]
				elif strand == '-':
					refpositions.reverse()
					ribosomeposition = refpositions[12]
			
			#Sometimes (rarely) ribosome position is None
			if ribosomeposition == None:
				continue

			consideredreadcounter +=1
			if chrm not in ribopositions:
				ribopositions[chrm] = {}
			if strand not in ribopositions[chrm]:
				ribopositions[chrm][strand] = []

			ribopositions[chrm][strand].append(ribosomeposition)

	print 'Found ribosome positions for {0} of {1} reads ({2}%).'.format(consideredreadcounter, readcounter, round((consideredreadcounter / float(readcounter)), 5)*100)

	return ribopositions

#########################
#########################
#Get longest CDS of gene#
#########################
#########################
def getCDScoords(gff):
	allCDScoords = {} #{chrm : {strand : [[list of CDS positions for gene 1], [list of CDS positions of gene2]]}}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	geneswithcodingtranscript = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1
		if genecount % 10000 == 0:
			print 'Gene {0}...'.format(genecount)
		geneboundaries[str(gene.id).replace('gene:', '')] = [gene.start, gene.end]

		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				geneswithcodingtranscript +=1
				break
		if not codingtranscript:
			continue


		CDSlengths = {} #{transcriptID : combined_length_of_coding_exons}
		CDScoords = {} #{transcriptID : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
		genename = str(gene.id)
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start'):
			transcriptID = str(transcript.id)
			CDScoords[transcriptID] = []
			CDSlength = 0
			for codingexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
				CDScoords[transcriptID].append([codingexon.start, codingexon.end])
				exonlength = codingexon.end - codingexon.start
				CDSlength += exonlength
			CDSlengths[transcriptID] = CDSlength

		longestcds = max(CDSlengths.iterkeys(), key = (lambda key: CDSlengths[key]))
		for transcript in CDScoords:
			if transcript == longestcds:
				cdsnt = []
				if chrm not in allCDScoords:
					allCDScoords[chrm] = {}
				if strand not in allCDScoords[chrm]:
					allCDScoords[chrm][strand] = []
				if strand == '+':
					for exon in CDScoords[transcript]:
						cdsnt += range(exon[0], exon[1] + 1)
				elif strand == '-':
					for exon in reversed(CDScoords[transcript]):
						cdsnt += list(reversed(range(exon[0], exon[1] + 1)))

				allCDScoords[chrm][strand].append(cdsnt)

	#os.remove(db_fn)
	#pickle.dump(allCDScoords, open('cds.p', 'wb'))
	return allCDScoords



def getlongest3putrcoords(gff):
	all3pUTRcoords = defaultdict(list) # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	longest3pUTRcoords = {} # {chrm : {strand : [[list of UTR positions for gene 1], [list of UTR positions of gene2]]}}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	no3pUTRcount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1
		if genecount % 5000 == 0:
			print 'Gene {0}...'.format(genecount)
		geneboundaries[str(gene.id).replace('gene:', '')] = [gene.start, gene.end]

		genename = str(gene.id) + '_' + str(gene.chrom) + '_' + gene.strand
		strand = gene.strand

		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				break
		if not codingtranscript:
			continue
		
		all3pUTRcoords[genename] #initiate entry as list (defaultdict)
		for transcript in db.children(gene, featuretype = 'transcript'):
			#If there's no CDS in this transcript, skip it
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				continue

			exoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
			CDScoords = []
			UTRcoords = [] #[UTRstart, UTRstop]
			for exon in db.children(transcript, featuretype = 'exon'):
				exoncoords.append([exon.start, exon.end])
			for CDSexon in db.children(transcript, featuretype = 'CDS'):
				CDScoords.append([CDSexon.start, CDSexon.end])

			#3' UTR start is directly after CDS end
			if transcript.strand == '+':
				CDSend = max(CDScoords, key = itemgetter(1))[1]
				if CDSend == transcript.end: #the transcript ends right where the CDS ends
					continue
				UTR3start = CDSend + 1
				UTRcoords = [UTR3start, transcript.end]
			elif transcript.strand == '-':
				CDSend = min(CDScoords, key = itemgetter(0))[0]
				if CDSend == transcript.start: #the transcript ends right where the CDS ends
					continue
				UTR3start = CDSend - 1
				UTRcoords = [transcript.start, UTR3start]

			
			###Check to see if the UTR is fully contained within the coordinates of one exon
			singleexonUTR = False
			for exoncoord in exoncoords:
				exonstart = exoncoord[0]
				exonend = exoncoord[1]
				if exonstart <= UTRcoords[0] and exonend >= UTRcoords[1]:
					singleexonUTR = True
					all3pUTRcoords[genename].append([[UTRcoords[0], UTRcoords[1]]])

			if singleexonUTR == False:
				#Get all positions that are both exonic and in the 3' UTR
				overlappingbp = [] #sorted exonic positions in UTR
				UTR3range = range(UTRcoords[0], UTRcoords[1] + 1)
				for exoncoord in exoncoords:
					exonrange = range(exoncoord[0], exoncoord[1] + 1)
					overlap = set(UTR3range).intersection(exonrange)
					for nt in sorted(list(overlap)):
						overlappingbp.append(nt)

				#Now get breaks in consecutive exonic positions
				#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
				UTRexoncoords = []
				for k, g in groupby(enumerate(overlappingbp), lambda (index, item): index-item):
					exonbp = map(itemgetter(1), g)
					if len(exonbp) > 1:
						UTRexoncoords.append([exonbp[0], exonbp[-1]])
				##############ADD FLIP OF EXON CONNECTIVITY FOR - STRAND HERE??????
				if strand == '-':
					UTRexoncoords.reverse()
				all3pUTRcoords[genename].append(UTRexoncoords)

	#Now get longest 3' UTR for each gene
	for gene in all3pUTRcoords:
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]
		if len(all3pUTRcoords[gene]) == 0: #if there were no 3' UTRs found for this gene
			continue
		UTRlengths = [] #list of lengths of UTRs for this gene, in order of appearance in all3pUTRcoords
		for UTR in all3pUTRcoords[gene]:
			#If this UTR has no coords, skip it
			if len(UTR) == 0:
				print gene, all3pUTRcoords[gene]
				continue
			if type(UTR[0]) == int:
				#This must be a single exon UTR
				UTRlength = UTR[1] - UTR[0]
				if UTRlength < 0:
					print 'Error! Negative UTR length for {0}.'.format(gene)
				UTRlengths.append(UTRlength)
			elif type(UTR[0]) == list:
				#This must be a multi-exon UTR
				UTRlength = 0
				for exon in UTR:
					exonlength = exon[1] - exon[0]
					if exonlength < 0:
						print 'Error! Negative UTR exon length for {0}'.format(gene)
					UTRlength += exonlength
				UTRlengths.append(UTRlength)

		indexofmax = max((length, index) for index, length in enumerate(UTRlengths))[1]
		#longest3pUTRcoords[gene] = all3pUTRcoords[gene][indexofmax]
		if chrm not in longest3pUTRcoords:
			longest3pUTRcoords[chrm] = {}
		if strand not in longest3pUTRcoords[chrm]:
			longest3pUTRcoords[chrm][strand] = []

		coords = []
		if strand == '+':
			for exon in all3pUTRcoords[gene][indexofmax]:
				coords += range(exon[0], exon[1] + 1)
				longest3pUTRcoords[chrm][strand].append(sorted(coords))
		elif strand == '-':
			for exon in all3pUTRcoords[gene][indexofmax]: #Exon order has already been reversed above
				coords += list(reversed(range(exon[0], exon[1] + 1)))
				longest3pUTRcoords[chrm][strand].append(list(reversed(sorted(coords))))



	#pickle.dump(longest3pUTRcoords, open('utr3.p', 'wb'))
	return longest3pUTRcoords

#################################
#Get longest 5' UTRs
#################################
def getlongest5putrcoords(gff):
	all5pUTRcoords = defaultdict(list) # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	longest5pUTRcoords = {} # {chrm : {strand : [[list of UTR positions for gene 1], [list of UTR positions of gene2]]}}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	no3pUTRcount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1
		if genecount % 5000 == 0:
			print 'Gene {0}...'.format(genecount)
		geneboundaries[str(gene.id).replace('gene:', '')] = [gene.start, gene.end]

		genename = str(gene.id) + '_' + str(gene.chrom) + '_' + gene.strand
		strand = gene.strand

		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				break
		if not codingtranscript:
			continue
		
		all5pUTRcoords[genename] #initiate entry as list (defaultdict)
		for transcript in db.children(gene, featuretype = 'transcript'):
			#If there's no CDS in this transcript, skip it
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				continue

			exoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
			CDScoords = []
			UTRcoords = [] #[UTRstart, UTRstop]
			for exon in db.children(transcript, featuretype = 'exon'):
				exoncoords.append([exon.start, exon.end])
			for CDSexon in db.children(transcript, featuretype = 'CDS'):
				CDScoords.append([CDSexon.start, CDSexon.end])

			#5' UTR end is directly before CDS start
			if transcript.strand == '+':
				CDSstart = min(CDScoords, key = itemgetter(0))[0]
				if CDSstart == transcript.start: #the transcript starts right where the CDS starts
					continue
				UTR5end = CDSstart - 1
				UTRcoords = [transcript.start, UTR5end]
			elif transcript.strand == '-':
				CDSstart = max(CDScoords, key = itemgetter(1))[1]
				if CDSstart == transcript.end: #the transcript starts right where the CDS starts
					continue
				UTR5end = CDSstart + 1
				UTRcoords = [UTR5end, transcript.end]

			
			###Check to see if the UTR is fully contained within the coordinates of one exon
			singleexonUTR = False
			for exoncoord in exoncoords:
				exonstart = exoncoord[0]
				exonend = exoncoord[1]
				if exonstart <= UTRcoords[0] and exonend >= UTRcoords[1]:
					singleexonUTR = True
					all5pUTRcoords[genename].append([[UTRcoords[0], UTRcoords[1]]])

			if singleexonUTR == False:
				#Get all positions that are both exonic and in the 3' UTR
				overlappingbp = [] #sorted exonic positions in UTR
				UTR5range = range(UTRcoords[0], UTRcoords[1] + 1)
				for exoncoord in exoncoords:
					exonrange = range(exoncoord[0], exoncoord[1] + 1)
					overlap = set(UTR5range).intersection(exonrange)
					for nt in sorted(list(overlap)):
						overlappingbp.append(nt)

				#Now get breaks in consecutive exonic positions
				#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
				UTRexoncoords = []
				for k, g in groupby(enumerate(overlappingbp), lambda (index, item): index-item):
					exonbp = map(itemgetter(1), g)
					if len(exonbp) > 1:
						UTRexoncoords.append([exonbp[0], exonbp[-1]])
				if strand == '-':
					UTRexoncoords.reverse()
				all5pUTRcoords[genename].append(UTRexoncoords)

	#Now get longest 5' UTR for each gene
	for gene in all5pUTRcoords:
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]
		if len(all5pUTRcoords[gene]) == 0: #if there were no 3' UTRs found for this gene
			continue
		UTRlengths = [] #list of lengths of UTRs for this gene, in order of appearance in all3pUTRcoords
		for UTR in all5pUTRcoords[gene]:
			if type(UTR[0]) == int:
				#This must be a single exon UTR
				UTRlength = UTR[1] - UTR[0]
				if UTRlength < 0:
					print 'Error! Negative UTR length for {0}.'.format(gene)
				UTRlengths.append(UTRlength)
			elif type(UTR[0]) == list:
				#This must be a multi-exon UTR
				UTRlength = 0
				for exon in UTR:
					exonlength = exon[1] - exon[0]
					if exonlength < 0:
						print 'Error! Negative UTR exon length for {0}'.format(gene)
					UTRlength += exonlength
				UTRlengths.append(UTRlength)

		indexofmax = max((length, index) for index, length in enumerate(UTRlengths))[1]
		coords = []
		#longest5pUTRcoords[gene] = all5pUTRcoords[gene][indexofmax]
		if chrm not in longest5pUTRcoords:
			longest5pUTRcoords[chrm] = {}
		if strand not in longest5pUTRcoords[chrm]:
			longest5pUTRcoords[chrm][strand] = []

		if strand == '+':
			for exon in all5pUTRcoords[gene][indexofmax]:
				coords += range(exon[0], exon[1] + 1)
				longest5pUTRcoords[chrm][strand].append(sorted(coords))
		elif strand == '-':
			for exon in all5pUTRcoords[gene][indexofmax]: #exon order was already reversed above
				coords += list(reversed(range(exon[0], exon[1] + 1)))
				longest5pUTRcoords[chrm][strand].append(list(reversed(sorted(coords))))


	#pickle.dump(longest5pUTRcoords, open('utr5.p', 'wb'))
	return longest5pUTRcoords


def binreads(ribopositions, UTR5coords, CDScoords, UTR3coords, numberofbins):
	#ribopositions = {} #{chrm : {strand : [positions]}}
	coorddicts = [UTR5coords, CDScoords, UTR3coords]
	bins = {} #{region : {binnumber : [readcount, maxpossiblereadcount]}}
	regionlengths = {} #{region : {chrm : {strand : {nt : [list of nts also in this region]}}}
	binassignments = {} #{region : {chrm : {strand : {position : bin}}}}
	regions = ['UTR5', 'CDS', 'UTR3']

	print 'Making length dict...'
	#Make length dict
	#How many bins exist for each CDS or UTR?
	#Not every bin may exist in every CDS/UTR due to rounding

	#For every position, have it point to a list of all other positions that are in the same CDS/UTR
	for idx, region in enumerate(regions):
		coorddict = coorddicts[idx]
		if region not in regionlengths:
			regionlengths[region] = {}
		for chrm in coorddict:
			if chrm not in regionlengths[region]:
				regionlengths[region][chrm] = {}
			for strand in coorddict[chrm]:
				if strand not in regionlengths[region][chrm]:
					regionlengths[region][chrm][strand] = {}
				for cds in coorddict[chrm][strand]:
					for nt in cds:
						regionlengths[region][chrm][strand][nt] = cds


	#Make bin assignmentdict
	#Assign each observed coordinate (from longest CDS, UTR, etc.) to a bin
	print 'Making bin assignments...'
	for indx, region in enumerate(regions):
		coorddict = coorddicts[indx]
		if region not in binassignments:
			binassignments[region] = {}
		for chrm in coorddict:
			if chrm not in binassignments[region]:
				binassignments[region][chrm] = {}
			for strand in coorddict[chrm]:
				if strand not in binassignments[region][chrm]:
					binassignments[region][chrm][strand] = {}
				for cds in coorddict[chrm][strand]:
					for idx, cdsnt in enumerate(cds):
						binnumber = math.floor((idx / float(len(cds))) * numberofbins)
						binassignments[region][chrm][strand][cdsnt] = binnumber

	
	#populate bins
	for region in regions:
		bins[region] = {}
		for binnumber in range(numberofbins):
			bins[region][binnumber] = [0, 0] #[readcount, maxpossiblereadcount]

	
	#Go through each read. Look to see if it was in a longest CDS, UTR, etc.
	#If not, just pass.
	#If it is, its bin assignment will be in binassignments
	#You can then add one "possible" bin assignment to every bin that would be hit in that region
	readcounter = 0
	binnedreadcounter = 0
	for chrm in ribopositions:
		for strand in ribopositions[chrm]:
			for riboposition in ribopositions[chrm][strand]:
				readcounter +=1
				if readcounter % 100000 == 0:
					print readcounter
				for region in regions:
					try:
						binnumber = binassignments[region][chrm][strand][riboposition]
						bins[region][binnumber][0] +=1
						#Get all possible binnumbers in this cds
						binnedreadcounter +=1
						for idx, othernt in enumerate(regionlengths[region][chrm][strand][riboposition]):
							possiblebinnumber = math.floor((idx / float(len(regionlengths[region][chrm][strand][riboposition]))) * numberofbins)
							bins[region][possiblebinnumber][1] +=1

					except KeyError:
						pass

	#for region in bins:
		#for binnumber in bins[region]:
			#print region, binnumber, bins[region][binnumber], bins[region][binnumber][0] / (float(bins[region][binnumber][1] + 1)), binnedreadcounter

	return bins, binnedreadcounter


#Instead of binning reads into bins, alternatively, just align them to UTR/CDS boundaries, up to 200 nt away from this boundary.
#nt 0 is the first nt of the region.
#nt -1 is the last nt of the region

def alignreads(ribopositions, UTR5coords, CDScoords, UTR3coords):
	#CDScoords = {} #{chrm : {strand : [[list of CDS positions for gene 1], [list of CDS positions of gene2]]}}
	#ribopositions = {} #{chrm : {strand : [positions]}}
	aligndict = {} #{region : {alignment : readcount}}
	positiondict = {} #{chrm : {strand : {position : assignment}}} #e.g. assignment = 'UTR3_-150'
	coorddicts = [UTR5coords, UTR3coords, CDScoords]
	regions = ['UTR5', 'UTR3', 'CDS']
	assignedreadcounter = 0 #number of reads that were assigned a -200:200 position

	
	#Make positiondict
	for ind, coorddict in enumerate(coorddicts):
		regionlabel = regions[ind]
		for chrm in coorddict:
			for strand in coorddict[chrm]:
				for region in coorddict[chrm][strand]:
					if len(region) < 400:
						continue
					for i in range(200):
						#Some regions will be smaller than 200 and will cause a problem
						try:
							nt = region[i]
						except IndexError:
							continue
						label = regionlabel + '_' + str(i)
						if chrm not in positiondict:
							positiondict[chrm] = {}
						if strand not in positiondict[chrm]:
							positiondict[chrm][strand] = {}
						positiondict[chrm][strand][nt] = label
					for i in range(-200, 0):
						try:
							nt = region[i]
						except IndexError:
							continue
						label = regionlabel + '_' + str(i)
						positiondict[chrm][strand][nt] = label


	#Populate aligndict
	for region in regions:
		aligndict[region] = {}
	for i in range(200):
		for region in regions:
			aligndict[region][i] = 0
	for i in range(-200, 0):
		for region in regions:
			aligndict[region][i] = 0

	#Assign ribopositions to alignments
	for chrm in ribopositions:
		for strand in ribopositions[chrm]:
			for riboposition in ribopositions[chrm][strand]:
				#Not every riboposition is going to have an entry in positiondict
				#because not every riboposition is within 200 nt of a boundary
				try:
					label = positiondict[chrm][strand][riboposition]
					region = label.split('_')[0]
					position = int(label.split('_')[1])
					aligndict[region][position] +=1
					assignedreadcounter +=1
				except KeyError:
					pass

	return aligndict, assignedreadcounter



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF annotation of genome. Used for defining CDS and UTR regions.')
	parser.add_argument('--bamdirectory', type = str, help = 'Directory containing bam files to examine.')
	parser.add_argument('--mode', type = str, choices = ['bins', 'align'], help = 'Make a bin-based metagene or align reads to UTR/CDS boundaries?')
	parser.add_argument('--numberofbins', type = int, help = 'Number of bins for metagene. Needed if mode = bins')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	
	if args.mode == 'bins':
		with open(args.outfile, 'w') as outfh:
			outfh.write(('\t').join(['sample', 'region', 'binnumber', 'readcount', 'maxpossiblereadcount', 'fraction', 'binnedreadcounter', 'normalizedfraction']) + '\n')

		CDScoords = getCDScoords(args.gff)
		UTR3coords = getlongest3putrcoords(args.gff)
		UTR5coords = getlongest5putrcoords(args.gff)

		for bamfile in os.listdir(args.bamdirectory):
			if bamfile.endswith('.bam'):
				sample = bamfile.split('.')[0]
				print sample
				ribopositions = getribopositions(bamfile)
				bins, binnedreadcounter = binreads(ribopositions, UTR5coords, CDScoords, UTR3coords, args.numberofbins)
				with open(args.outfile, 'a') as outfh:
					for region in bins:
						for binnumber in bins[region]:
							assignedreads = bins[region][binnumber][0]
							maxpossiblereadcount = bins[region][binnumber][1] + 1 #pseudocount
							fraction = assignedreads / float(maxpossiblereadcount)
							#Normalize to number of reads
							normalizedfraction = (fraction / float(binnedreadcounter)) * 1000000
							outfh.write(('\t').join([sample, region, str(binnumber), str(assignedreads), str(maxpossiblereadcount), str(fraction), str(binnedreadcounter), str(normalizedfraction)]) + '\n')

	elif args.mode == 'align':
		with open(args.outfile, 'w') as outfh:
			outfh.write(('\t').join(['sample', 'region', 'alignment', 'readcount', 'samplereadcount', 'normalizedreadcount']) + '\n')

		CDScoords = getCDScoords(args.gff)
		UTR3coords = getlongest3putrcoords(args.gff)
		UTR5coords = getlongest5putrcoords(args.gff)

		for bamfile in os.listdir(args.bamdirectory):
			if bamfile.endswith('.bam'):
				sample = bamfile.split('.')[0] + '.' + bamfile.split('.')[1]
				print sample
				ribopositions = getribopositions(bamfile)
				alignpositions, assignedreadcount = alignreads(ribopositions, UTR5coords, CDScoords, UTR3coords)
				with open(args.outfile, 'a') as outfh:
					for region in alignpositions:
						for alignposition in alignpositions[region]:
							readcount = alignpositions[region][alignposition]
							normalizedreadcount = (readcount / float(assignedreadcount)) * 1000000
							outfh.write(('\t').join([sample, region, str(alignposition), str(readcount), str(assignedreadcount), str(normalizedreadcount)]) + '\n')



