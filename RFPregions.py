#Where do the RFP reads map to? CDS? 3' UTR? 5' UTR? intron? intergenic?
#Each read is going to have one position that it "aligns" to.  This will not, therefore,
#be coverage.  The reason for doing it this way is because the place that the ribosome
#is "translating" is 15 nt offset from the beginning of the read.  Therefore, the beginning
#of the read could, for example, be in a 5' UTR, while the ribosome is really doing its business
#in the CDS.

import pysam
import gffutils
import os
from operator import itemgetter
from itertools import groupby
import argparse

#Given a bam, get the ribosome "positions" (15 nt offset from read start)
#Importantly, due to introns, this may not just be readstart + 15!
def getribopositions(bam):
	readcounter = 0
	consideredreadcounter = 0
	ribopositions = {} #{chrm : {strand : [positions]}}
	if 'tot' in bam:
		acceptedreadlengths = [44, 45, 46, 47, 48, 49, 50]
	elif 'rpf' in bam or 'RFP' in bam:
		acceptedreadlengths = [30, 31, 32, 33, 34, 35]
	elif 'rnaseq' in bam:
		acceptedreadlengths = [40]
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
			#Then, take the position that is right in the middle of that list as the "position" of the ribosome
			halfway = (read.query_length / 2) - 1 #make it zero based
			refpos = read.get_reference_positions(full_length = True) #these are 0-based since they are coming from a bam
			ribosomeposition = refpos[halfway]
			
			#Sometimes (rarely) ribosome position is None
			if ribosomeposition == None:
				continue

			consideredreadcounter +=1
			if chrm not in ribopositions:
				ribopositions[chrm] = {}
			if strand not in ribopositions[chrm]:
				ribopositions[chrm][strand] = []

			#Turn them into 1-based coords for easy matching with gff
			ribosomeposition = ribosomeposition + 1

			ribopositions[chrm][strand].append(ribosomeposition)

	print 'Found ribosome positions for {0} of {1} reads ({2}%).'.format(consideredreadcounter, readcounter, round((consideredreadcounter / float(readcounter)), 5)*100)

	return ribopositions


#Given a gff, assign every genic position to either CDS, 3' UTR, 5' UTR, intron, or intergenic.
#There are obviously going to be some overlaps. Something that is 3' UTR in one transcript may be
#CDS in another.  Therefore, we have to have hierarchies.  If something is CDS in one transcript,
#it will be annotated as such no matter what it is in any other transcript.
#The hierarchy is CDS > 3' UTR > 5' UTR > noncodingexon > intron > intergenic.

#First, get all CDS coords
def getCDScoords(gff):
	allCDScoords = {} #{chrm : {strand : [list of CDS positions]}}
	genecount = 0
	geneswithcodingtranscript = 0
	


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
		chrm = str(gene.chrom)
		strand = str(gene.strand)
		if chrm not in allCDScoords:
			allCDScoords[chrm] = {}
		if strand not in allCDScoords[chrm]:
			allCDScoords[chrm][strand] = []

		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				geneswithcodingtranscript +=1
				break
		if not codingtranscript:
			continue

		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			transcriptid = str(transcript.id)
			for codingexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
				allCDScoords[chrm][strand] += range(codingexon.start, codingexon.end + 1)

	#Remove duplicate positions
	for chrm in allCDScoords:
		for strand in allCDScoords[chrm]:
			positions = allCDScoords[chrm][strand]
			uniquepositions = sorted(list(set(positions)))
			allCDScoords[chrm][strand] = uniquepositions

	return allCDScoords

#Now get any coord that is in a 3' UTR in any transcript
def getUTR3coords(gff):
	allUTR3coords = {} #{chrm : {strand : [list of 3' UTR coords]}}
	genecount = 0

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
		chrm = str(gene.chrom)
		strand = str(gene.strand)

		if chrm not in allUTR3coords:
			allUTR3coords[chrm] = {}
		if strand not in allUTR3coords[chrm]:
			allUTR3coords[chrm][strand] = []

		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			#If theres no CDS in this transcript, skip it
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
					UTRexoncoords = [[UTRcoords[0], UTRcoords[1]]]

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
				UTRexoncoords = [] #[[UTRexon1start, UTRexon1end], [UTRexon2start, UTRexon2end]]
				for k, g in groupby(enumerate(overlappingbp), lambda (index, item): index-item):
					exonbp = map(itemgetter(1), g)
					if len(exonbp) > 1:
						UTRexoncoords.append([exonbp[0], exonbp[-1]])
				##############ADD FLIP OF EXON CONNECTIVITY FOR - STRAND HERE??????
				if strand == '-':
					UTRexoncoords.reverse()

			for UTRexon in UTRexoncoords:
				allUTR3coords[chrm][strand] += range(UTRexon[0], UTRexon[1] + 1)

	#Remove duplicate positions
	for chrm in allUTR3coords:
		for strand in allUTR3coords[chrm]:
			positions = allUTR3coords[chrm][strand]
			uniquepositions = sorted(list(set(positions)))
			allUTR3coords[chrm][strand] = uniquepositions

	return allUTR3coords

def getUTR5coords(gff):
	allUTR5coords = {} #{chrm : {strand : [list of 3' UTR coords]}}
	genecount = 0

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
		chrm = str(gene.chrom)
		strand = str(gene.strand)

		if chrm not in allUTR5coords:
			allUTR5coords[chrm] = {}
		if strand not in allUTR5coords[chrm]:
			allUTR5coords[chrm][strand] = []
				
		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				break
		if not codingtranscript:
			continue

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
					UTRexoncoords = [[UTRcoords[0], UTRcoords[1]]]

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

			for UTRexon in UTRexoncoords:
				allUTR5coords[chrm][strand] += range(UTRexon[0], UTRexon[1] + 1)

	#Remove duplicate positions
	for chrm in allUTR5coords:
		for strand in allUTR5coords[chrm]:
			positions = allUTR5coords[chrm][strand]
			uniquepositions = sorted(list(set(positions)))
			allUTR5coords[chrm][strand] = uniquepositions

	return allUTR5coords

#Get the exonic coords of lncrnas
def getnoncodingcoords(gff):
	allnccoords = {} #{chrm : {strand : [list of 3' UTR coords]}}
	genecount = 0

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
		chrm = str(gene.chrom)
		strand = str(gene.strand)

		if chrm not in allnccoords:
			allnccoords[chrm] = {}
		if strand not in allnccoords[chrm]:
			allnccoords[chrm][strand] = []

		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			#We want transcripts that have no annotated CDS
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				for exon in db.children(transcript, featuretype = 'exon', level = 1):
					allnccoords[chrm][strand] += range(exon.start, exon.end + 1)

	for chrm in allnccoords:
		for strand in allnccoords[chrm]:
			positions = allnccoords[chrm][strand]
			uniquepositions = sorted(list(set(positions)))
			allnccoords[chrm][strand] = uniquepositions

	return allnccoords

def getintroncoords(gff):
	allintroncoords = {} #{chrm : {strand : [list of intron coords]}}

	genecount = 0

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
		chrm = str(gene.chrom)
		strand = str(gene.strand)

		if chrm not in allintroncoords:
			allintroncoords[chrm] = {}
		if strand not in allintroncoords[chrm]:
			allintroncoords[chrm][strand] = []

		genecoords = range(gene.start, gene.end + 1)
		exoncoords = []
		for exon in db.children(gene, featuretype = 'exon'):
			exoncoords += range(exon.start, exon.end + 1)

		introncoords = set(genecoords) - set(exoncoords)
		introncoords = sorted(list(set(introncoords)))

		allintroncoords[chrm][strand] += introncoords

	return allintroncoords


def assignpositions(gff, introncoords, noncodingcoords, UTR5coords, UTR3coords, CDScoords):
	assignments = {} #{chrm : {strand : {position : <'CDS' or 'UTR3' or 'UTR5' or 'intron' or 'intergenic'>}}}
	hierarchy = {'CDS' : 6, 'UTR3' : 5, 'UTR5' : 4, 'noncodingexon' : 3, 'intron' : 2, 'intergenic' : 1}
	chrms = [] #all chromosomes in gff
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		chrm = str(gene.chrom)
		if chrm not in chrms:
			chrms.append(chrm)

	#Populate assignments
	for chrm in chrms:
		assignments[chrm] = {}
		for strand in ['+', '-']:
			assignments[chrm][strand] = {}

	#Add intronic coords
	for chrm in introncoords:
		for strand in introncoords[chrm]:
			for position in introncoords[chrm][strand]:
				assignments[chrm][strand][position] = 'intron'

	#Add noncoding coords
	for chrm in noncodingcoords:
		for strand in noncodingcoords[chrm]:
			for position in noncodingcoords[chrm][strand]:
				assignments[chrm][strand][position] = 'noncoding'

	#Add UTR5 coords
	for chrm in UTR5coords:
		for strand in UTR5coords[chrm]:
			for position in UTR5coords[chrm][strand]:
				assignments[chrm][strand][position] = '5UTR'

	#Add UTR3 coords
	for chrm in UTR3coords:
		for strand in UTR3coords[chrm]:
			for position in UTR3coords[chrm][strand]:
				assignments[chrm][strand][position] = '3UTR'

	#Add CDS coords
	for chrm in CDScoords:
		for strand in CDScoords[chrm]:
			for position in CDScoords[chrm][strand]:
				assignments[chrm][strand][position] = 'CDS'

	return assignments

def assignribopositions(ribopositions, assignments):
	assignmentdict = {'CDS' : 0, '3UTR' : 0, '5UTR' : 0, 'noncoding' : 0, 'intron' : 0, 'intergenic' : 0}
	for chrm in ribopositions:
		for strand in ribopositions[chrm]:
			for position in ribopositions[chrm][strand]:
				try:
					assignment = assignments[chrm][strand][position]
				except KeyError:
					assignment = 'intergenic'

				assignmentdict[assignment] +=1

	return assignmentdict


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF of genome annotation.')
	#parser.add_argument('--bam', type = str, help = 'BAM file of rpf alignments.')
	parser.add_argument('--bamdirectory', type = str, help = 'Bam directory.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	
	print 'Finding CDS coords...'
	cdscoords = getCDScoords(args.gff)
	print 'Finding 3UTR coords...'
	UTR3coords = getUTR3coords(args.gff)
	print 'Finding 5UTR coords...'
	UTR5coords = getUTR5coords(args.gff)
	print 'Finding noncoding exon coords...'
	noncodingcoords = getnoncodingcoords(args.gff)
	print 'Finding intronic coords...'
	introncoords = getintroncoords(args.gff)
	

	assignments = assignpositions(args.gff, introncoords, noncodingcoords, UTR5coords, UTR3coords, cdscoords)
	with open(args.outfile, 'w') as outfh:
		outfh.write(('\t').join(['sample', 'region', 'fraction']) + '\n')
		for bamfile in os.listdir(args.bamdirectory):
			print bamfile
			if bamfile.endswith('.bam'):
				ribopositions = getribopositions(bamfile)
				assignmentdict = assignribopositions(ribopositions, assignments)
				totalreads = float(sum(assignmentdict.values()))
				for region in assignmentdict:
					fraction = assignmentdict[region] / totalreads
					outfh.write(('\t').join([bamfile, region, str(fraction)]) + '\n')

