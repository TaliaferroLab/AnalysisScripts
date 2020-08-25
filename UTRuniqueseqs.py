import gffutils
import os
from collections import defaultdict
import sys
from itertools import groupby, chain
from operator import itemgetter
import pickle


def transcriptfilters(transcript, db):
	mrnaendpass = False
	transcriptsupportlevelpass = False

	#Are we confident in the 3' end of this mrnaendpass
	if 'tag' not in transcript.attributes or 'mRNA_end_NF' not in transcript.attributes['tag']:
		mrnaendpass = True

	try:
		tsl = transcript.attributes['transcript_support_level'][0]
	except KeyError: #if this transcript doesnt have a transcript support level
		tsl = None
	if tsl != 'NA':
		transcriptsupportlevelpass = True

	if mrnaendpass and transcriptsupportlevelpass:
		return True
	else:
		return False


def getall3pUTR(gff):
	all3pUTRcoords = defaultdict(list) # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	genecount = 0
	geneboundaries = {} # {ensid : [genestart, genestop]}
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
		for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'start'):
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

			#Does this transcript pass filters?
			passtranscriptfilters = transcriptfilters(transcript, db)
			if passtranscriptfilters == False:
				continue

			exoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
			CDScoords = []
			UTRcoords = [] #[UTRstart, UTRstop]
			for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
				exoncoords.append([exon.start, exon.end])
			for CDSexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
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

					#No duplicates
					if [UTRcoords[0], UTRcoords[1]] not in all3pUTRcoords[genename]:
						all3pUTRcoords[genename].append([UTRcoords[0], UTRcoords[1]])


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
				#if strand == '-':
					#UTRexoncoords.reverse()

				#No duplicates
				if UTRexoncoords not in all3pUTRcoords[genename]:
					all3pUTRcoords[genename].append(UTRexoncoords)

	#Remove empty genes
	#These are genes that have no UTRs assigned to them.
	all3pUTRcoords = {gene:all3pUTRcoords[gene] for gene in all3pUTRcoords if all3pUTRcoords[gene]}


	print 'Looked through {0} genes and retrived UTRs for {1} of them.'.format(genecount, len(all3pUTRcoords))
	return all3pUTRcoords

def orderUTRs(all3pUTRcoords):
	#all3pUTRcoords = {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	all3pUTRcoords_sorted = {} # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	#For each gene:
	#order UTRs by end (max end for + strand and min end for - strand)

	for gene in all3pUTRcoords:
		genename = gene.split('_')[0]
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]

		utrs = all3pUTRcoords[gene]
		#Skip genes that only have one UTR
		if len(utrs) == 1:
			continue

		#print genename, strand
		#print utrs

		#Get the endpoints for each UTR
		utrends = [] # list of utr end coords for this gene
		if strand == '+':
			for utr in utrs:
				#If this is a single exon UTR, this will NOT be a nested list
				if any(isinstance(i, list) for i in utr) == False:
					utrend = max(utr)
					utrends.append(utrend)
				#If this is not a single exon UTR, this will be a nested list
				#If this is a multi-exon UTR, it's slightly more complicated, but still the maximum value
				#(doesn't matter if it's start or stop because stop is always greater than start) is what we want
				elif any(isinstance(i, list) for i in utr) == True:
					coords = list(chain.from_iterable(utr))
					utrend = max(coords)
					utrends.append(utrend)

		elif strand == '-':
			for utr in utrs:
				#If this is a single exon UTR, this will NOT be a nested list
				if any(isinstance(i, list) for i in utr) == False:
					utrend = min(utr)
					utrends.append(utrend)
				#If this is not a single exon UTR, this will be a nested list
				#If this is a multi-exon UTR, it's slightly more complicated, but still the maximum value
				#(doesn't matter if it's start or stop because stop is always greater than start) is what we want
				elif any(isinstance(i, list) for i in utr) == True:
					coords = list(chain.from_iterable(utr))
					utrend = min(coords)
					utrends.append(utrend)

		#print utrends



		#Sort UTRs by their ends
		#This can be done by sorting UTRs according to the values in utrends, since the ordering of the two is the same.
		if strand == '+':
			utrs_ends = zip(utrs, utrends)
			utrs_ends.sort()
			#print utrs_ends
			utrs_sorted = [u for u, e in utrs_ends]
		elif strand == '-':
			utrs_ends = zip(utrs, utrends)
			utrs_ends.sort(reverse = True)
			#print utrs_ends
			utrs_sorted = [u for u, e in utrs_ends]

		all3pUTRcoords_sorted[gene] = utrs_sorted

	return all3pUTRcoords_sorted



#all3pUTRcoords = getall3pUTR(sys.argv[1])
#outputfh = open('data.pkl', 'wb')
#pickle.dump(all3pUTRcoords, outputfh)
#outputfh.close()
pkl = open('data.pkl', 'rb')
all3pUTRcoords = pickle.load(pkl)
pkl.close()
all3pUTRcoords_sorted = orderUTRs(all3pUTRcoords)
print all3pUTRcoords_sorted


