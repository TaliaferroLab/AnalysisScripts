import gffutils
import os
import sys
import pysam
import subprocess
import argparse
import cPickle as pickle
import pybedtools
from Bio import SeqIO
import gzip
from itertools import product
from random import randrange
import scipy.stats as stats
import statsmodels.stats.multitest as mt
from collections import OrderedDict
from math import log
from random import sample

#bedtools 2.24.0
#pysam 0.15.0
#pybedtools 0.7.9


#Given an annotation as gff, get the longest CDS for every gene.
#Then get lists of nucleotide positions for every CDS for the purpose of assigning a frame to every nt in a "longest CDS".
def getCDScoords(gff):
	allCDScoords = {} #{ENSGENE_chrm_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	allCDSnt = {} #{chrm : {strand : [[list of nucleotides of CDS1], [list of nucleotides of CDS2]]}}
	framedict = {} #{chrm : {strand : {position (1-based) : frame}}}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	geneswithcodingtranscript = 0
	e2sdict = {} #{ENSGene : shortname}

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
				allCDScoords[genename + '_' + chrm + '_' + strand] = CDScoords[transcript]

	#os.remove(db_fn)
	#Now reorganize allCDScoords into allCDSnt
	#When this happens, do we want to remove the first 50 (or first 51 to keep frame easy) to minimize cycloheximide artefacts?
	for gene in allCDScoords:
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]
		if chrm not in allCDSnt:
			allCDSnt[chrm] = {}
		if strand not in allCDSnt[chrm]:
			allCDSnt[chrm][strand] = []

		cdsnt = []
		for exon in allCDScoords[gene]:
			exonnt = range(exon[0], exon[1] + 1)
			cdsnt += exonnt

		#Flip the coords for minus strand
		if strand == '-':
			cdsnt = cdsnt[::-1]

		#Remove the first 60 nt of the cds
		cdsnt = cdsnt[60:] 

		#Some CDS coords (about 2%) do not have a multiple of 3 length. They tend to be for not well-annotated genes (e.g. Gm20946)
		if len(cdsnt) %3 == 0:
			allCDSnt[chrm][strand].append(cdsnt)

	print 'Looked through {0} genes. {1} of them had non-NMD coding transcripts. Found longest CDS sequences for {2} of them.'.format(genecount, geneswithcodingtranscript, len(allCDScoords))

	cdsnt = 0
	nonframe = 0
	for chrm in allCDSnt:
		for strand in allCDSnt[chrm]:
			for cds in allCDSnt[chrm][strand]:
				if len(cds) % 3 != 0:
					nonframe +=1
				for nt in cds:
					cdsnt +=1

	#Now populate framedict
	for chrm in allCDSnt:
		if chrm not in framedict:
			framedict[chrm] = {}
		for strand in allCDSnt[chrm]:
			if strand not in framedict[chrm]:
				framedict[chrm][strand] = {}
			for CDS in allCDSnt[chrm][strand]:
				for ind, position in enumerate(CDS):
					frame = ind % 3
					framedict[chrm][strand][position] = frame

	return framedict

#Make a bed file where each line in the start/stop of the cds for the longest cds transcript of each gene
#we aren't caring about introns in the output bed
def makecdsbedandfasta(gff, genomefasta):
	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	print 'Done indexing!'

	genes = db.features_of_type('gene')
	genecount = 0
	with open('longestcds.bed', 'w') as bedoutfh, open('longestcds.fa', 'w') as fastaoutfh:
		for gene in genes:
			genecount +=1
			if genecount % 10000 == 0:
				print 'Gene {0}...'.format(genecount)

			#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
			codingtranscript = False
			for transcript in db.children(gene, featuretype = 'transcript', level = 1):
				if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
					codingtranscript = True
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
			longestcdsstart = CDScoords[longestcds][0][0]
			longestcdsend = CDScoords[longestcds][-1][1]

			#Write sequence
			longestcdsexons = CDScoords[longestcds]
			cdsseq = ''
			for exon in longestcdsexons:
				start = exon[0]
				end = exon[1]
				if strand == '+':
					exonseq = seq_dict[chrm].seq[start-1:end].upper()
					cdsseq += exonseq
				elif strand == '-':
					exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
					newseq = exonseq + cdsseq
					cdsseq = newseq

			#Only want longest CDSes that are divisible by 3
			if len(cdsseq) %3 != 0:
				#print 'ERROR: CDS of gene {0} is not divisible by 3.'.format(gene.id.split('.')[0])
				continue

			bedoutfh.write(('\t').join([chrm, str(longestcdsstart), str(longestcdsend), gene.id.split('.')[0], '1000', strand]) + '\n')
			fastaoutfh.write('>' + gene.id.split('.')[0] + '\n' + str(cdsseq) + '\n')


#To compare P site densities, we need to see what the expected background is.
#This is the codon frequency for every gene for every read that we consider.
#Every time a gene has a read that is mapped to it, one "iteration" of its background frequencies go into the background pool.
def getbackgroundcodonfreqs(longestcdsfasta):
	backgroundcodoncounts = {} #{geneid : {codon : count}}
	backgroundcodonfreqs = {} #{geneid : {codon : count / sum(all codon counts)}}
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]
	for record in SeqIO.parse(longestcdsfasta, 'fasta'):
		backgroundcodoncounts[record.id] = {}
		for codon in allcodons:
			backgroundcodoncounts[record.id][codon] = 0
		#Go through the seq, for every codon count how many times we see it
		seq = str(record.seq)
		for i in range(len(seq))[::3]: 
			codon = seq[i : i + 3]
			backgroundcodoncounts[record.id][codon] += 1

	#Turn counts into frequencies
	for gene in backgroundcodoncounts:
		backgroundcodonfreqs[gene] = {}
		totalcodoncounts = sum(backgroundcodoncounts[gene].values())
		for codon in backgroundcodoncounts[gene]:
			freq = backgroundcodoncounts[gene][codon] / float(totalcodoncounts)
			backgroundcodonfreqs[gene][codon] = freq


	return backgroundcodonfreqs, backgroundcodoncounts


#Given a read, get its P site codon.
def getcodon(read, positionframedict):
	if read.is_reverse:
		strand = '-'
	elif not read.is_reverse:
		strand = '+'
	chrm = read.reference_name
	#No mitochondrial ribosome stuff
	if chrm == 'chrM':
		return None, None
	if strand == '+':
		readstart = read.reference_start + 1
	elif strand == '-':
		readstart = read.reference_end

	#Offset readstart by 13 to get P site
	refpositions = read.get_reference_positions(full_length = True)
	#Using full_length = True will allow soft-clipped nts to have placeholders of None in this list
	#We still want them as placeholders so that the offset will be correct
	refpositions = [pos + 1 if type(pos) == int else pos for pos in refpositions] #make 1 based
	#if read.query_length == 30 or read.query_length == 31 or read.query_length == 32 or read.query_length == 33 or read.query_length == 34 or read.query_length == 35:
	if read.query_length == 31 or read.query_length == 32:
		if strand == '+':
			refpos = refpositions[13]
			codonstart = 13
		elif strand == '-':
			refpositions.reverse()
			refpos = refpositions[13]
			codonstart = 13

	elif read.query_length == 28 or read.query_length == 29 or read.query_length == 30:
		if strand == '+':
			refpos = refpositions[12]
			codonstart = 12
		elif strand == '-':
			refpositions.reverse()
			refpos = refpositions[12]
			codonstart = 12	

	#Not all readstarts are going to be in positionframedict.
	#Not all reads map to CDS.
	try:
		frame = positionframedict[chrm][strand][refpos]
		#If the P site for this read is in frame 0, get the codon in the x site
		if frame == 0:
			if strand == '+':
				readseq = read.query_sequence
			elif strand == '-':
				readseq = read.get_forward_sequence()
			codon = readseq[codonstart:codonstart + 3]
		elif frame != 0:
			codonstart, codon = None, None
	except KeyError:
		codonstart, codon = None, None

	return codonstart, codon

#Given a read, found out which gene it maps to.
def getgene(read, cdsbedtool):
	#Convert read to bed
	chrm = read.reference_name
	readstart = read.reference_start
	readend = read.reference_end
	readname = read.query_name
	if read.is_reverse:
		strand = '-'
	elif not read.is_reverse:
		strand = '+'
	bedline = [chrm, str(readstart), str(readend), readname, '1000', strand]
	with open('temp.bed', 'w') as outfh:
		outfh.write(('\t').join(bedline))
	readbedtool = pybedtools.BedTool('temp.bed')

	genematches = cdsbedtool.intersect(readbedtool, s = True)
	
	if len(genematches) == 1: #if this read only matches to one gene (or more precisely, longest CDS of one gene)
		gene = str(genematches).split('\t')[3]
	elif len(genematches) != 1:
		gene = None

	os.remove('temp.bed')
	return gene

#Make a dictionary of all CDS sequences
def makeseqdict():
	longestcdsdict = {} #{geneid : longestcdsseq}
	with open('longestcds.fa', 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			geneid = str(record.id)
			geneseq = str(record.seq)
			longestcdsdict[geneid] = geneseq

	return longestcdsdict


def getgenestringmatch(read, geneseq, psitecodonstart, psitecodon):
	#Given a read sequence and the sequence of the CDS that it maps to, figure out where within that gene that the alignment occurs
	#This allows you to get codons upstream and downstream of the p site
	#Psitecodonstart is the position in the *read* where the psite starts

	#For psite, codon index = 0
	#For asite, codon index = 1
	#For esite, codon index = -1
	surroundingcodons = {} #{codonindex : codon}

	if read.is_reverse:
		strand = '-'
	elif not read.is_reverse:
		strand = '+'
	
	#Get sequence of the read
	if strand == '+':
		readseq = read.query_sequence
	elif strand == '-':
		readseq = read.get_forward_sequence()

	#Get where in the gene CDS this read maps
	try:
		readindex = geneseq.index(readseq)
	except ValueError:
		#Couldn't find a perfect match for this read in the gene (could be sequencing error, etc.)
		return None

	psiteindex = readindex + psitecodonstart #position in the *gene* where the psite starts
	psitecodoningene = geneseq[psiteindex : psiteindex + 3] #sequence of the derived psite in the *gene* for checking purposes

	#print psitecodon, psitecodoningene #hopefully these two are the same

	#We are only doing this if the read was in frame, so we can go upstream 300 nt and assume that that is 100 codons
	#Window we are interested in is 300 nt upstream of psite start and 300 nt downstream of p site end
	windowstart = max(60, psiteindex - 300) #don't consider the first 60 nt of the CDS
	windowend = min(len(geneseq), psiteindex + 3 + 300)
	windowstartindex = int((psiteindex - windowstart) / 3)
	windowendindex = int((windowend - psiteindex - 3) / 3)

	windowseq = geneseq[windowstart : windowend]
	if len(windowseq) % 3 != 0: #not sure how this could be true
		print('WINDOWSEQ NOT DIVISIBLE BY 3')
		return None

	#Take the windowseq, divide it up into codons, and add it to surroundingcodons
	i = windowstart
	windowseqindex = 0
	while i < windowend:
		#codonindex is the position of the codon relative to the p site
		codonindex = (i - psiteindex) / 3
		currentcodon = windowseq[windowseqindex : windowseqindex + 3]
		surroundingcodons[codonindex] = currentcodon
		if codonindex == 0: #if this is p site
			with open('windowseqs.fa', 'a') as outfh:
				outfh.write('>' + str(read.query_name) + '_' + str(psiteindex) + '_' + psitecodon + '_' + currentcodon + '\n' + windowseq + '\n')
		i += 3
		windowseqindex += 3

	return surroundingcodons


#Given a set of gene coordinates (or CDS coords) and a bam, get the reads that map to each gene
#Also make a dictionary telling us which gene each read maps to
def getcountspergene_tabix(bam, cdsbed):
	genecounts = {} #{gene : [readIDs that map to it]}
	reads2gene = {} #{readid : gene}

	#Make tabix index of the bam, and to do that first we have to convert to sam
	print 'Converting bam to sam...'
	with open('temp.sam', 'w') as outfh:
		command = ['samtools', 'view', '-h', bam]
		subprocess.call(command, stdout = outfh)

	#Compress sam
	print 'Compressing sam...'
	command = ['bgzip', 'temp.sam']
	subprocess.call(command)

	#Create tabix index
	print 'Creating tabix index...'
	command = ['tabix', '-p', 'sam', 'temp.sam.gz']
	subprocess.call(command)


	tbx = pysam.TabixFile('temp.sam.gz')
	with open(cdsbed, 'r') as cdsbedfh:
		genecounter = 0
		for line in cdsbedfh:
			genecounter +=1
			if genecounter % 5000 == 0:
				print 'Gene {0}...'.format(genecounter)
			line = line.strip().split('\t')
			chrm = line[0]
			start = int(line[1])
			end = int(line[2])
			genename = line[3]
			genebed = (' ').join(line)

			genecounts[genename] = []
			#If there is a gene that maps to a chromosome not present in any read (chrY, for example)
			#this will cause an error.  Except that error.
			try:
				readmatches = tbx.fetch(chrm, start, end)
			except ValueError:
				continue
			for readmatch in readmatches:
				#print readmatch
				readid = str(readmatch).split('\t')[0]
				reads2gene[readid] = genename
				genecounts[genename].append(readid)

	with open('{0}genecounts.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(genecounts, outfh)

	with open('{0}reads2gene.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(reads2gene, outfh)

	os.remove('temp.sam.gz')
	os.remove('temp.sam.gz.tbi')

	return genecounts, reads2gene




#Write a new bam containing all reads for which we can derive a valid p site codon
#This means that the read maps to one gene, and it's P site is in frame.
def getvalidreads(bam, positionframedict, cdsbed):
	#First, narrow down the bam to the reads we want to work with
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]
	with pysam.AlignmentFile(bam, 'rb') as infh, pysam.AlignmentFile('{0}validreads.bam'.format(os.path.basename(bam)), 'wb', template = pysam.AlignmentFile(bam, 'rb')) as outfh:
		readcounter = 0
		filteredreadcounter = 0
		for read in infh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print 'Read {0}...'.format(readcounter)
			#Only consider reads of lengths 28-32
			#if read.query_length not in [28, 29, 30, 31, 32]:
			if read.query_length not in [28, 29, 30]:
				continue

			#For non-RFP reads, the first read of mate pair is on the wrong strand
			#Also, only consider unique mappers
			if read.is_read1 or read.get_tag('NH') > 1 or read.is_secondary:
				continue


			psitecodonstart, psitecodon = getcodon(read, positionframedict)
			if psitecodon in allcodons:
				outfh.write(read)
				filteredreadcounter +=1

	print '{0} reads in this bam. {1} passed filters.'.format(readcounter, filteredreadcounter)


def iteratereads(filteredbam, positionframedict, cdsbed, reads2gene, longestcdsdict):
	with open('windowseqs.fa', 'w') as outfh:
		pass

	codondict = {} #{codonindex : {codon : [list of readIDs that have this codon at this index]}}
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]
	for codonindex in range(-100, 101):
		codondict[codonindex] = {}
		for codon in allcodons:
			codondict[codonindex][codon] = []


	cdsbedtool = pybedtools.BedTool(cdsbed)
	readcounter = 0

	#First, narrow down the bam to the reads we want to work with
	#Shouldn't have to do this anymore because we will have prefiletered the bam with getvalidreads()
	with pysam.AlignmentFile(filteredbam, 'rb') as infh:

		for read in infh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print 'Read {0}...'.format(readcounter)

			psitecodonstart, psitecodon = getcodon(read, positionframedict)
			gene = reads2gene[str(read.query_name)]
			geneseq = longestcdsdict[gene]
			surroundingcodons = getgenestringmatch(read, geneseq, psitecodonstart, psitecodon) #{codonindex : codon}
			#We have been returned a dictionary of all codons at the -100 to +100 sites relative to the p site.
			if surroundingcodons:
				for codonindex in surroundingcodons:
					codon = surroundingcodons[codonindex]
					codondict[codonindex][codon].append(read.query_name)

	return codondict


#Compare the frequencies for a codon being in the P site to expected background across all genes
#These frequencies are weighted by the number of times a gene had a Psite read in it
#Each gene has its own bg freq (that's in backgroundcodonfreqs).
#These freqs are multiplied by the number of Psite reads in that gene (then at the end divided by the total number of P site reads to get back to freqs).
def comparefreqs(codondict, genecountdict, backgroundcodonfreqs, backgroundcodoncounts):
	#codondict = {} #{codonindex : {codon : [list of readIDs that have this codon at this index]}}
	#genecountdict = {} #{gene : [readIDs that map to it]}
	#backgroundcodonfreqs = {} #{geneid : {codon : freq in that gene}}
	#backgroundcodoncounts = {} #{geneid : {codon : count}}
	totalbackgroundcounts = {} #{codon : count} (this is weighted by genecount)
	totalbackgroundcodondict = {} #{codon: final expected bg freq (all of these should add up to 1)}
	codoncounts = {} #{codonindex : {codon : [number of counts in xsite, number of total counts (all codons), number of counts in background, number of total counts (all codons) in background]}}
	enrichments = {} #{codonindex :{codon : log2 enrichment over background}}
	pvalues = OrderedDict() #{codonindex : {codon : pvalue}}
	adjustedpvalues = OrderedDict() #{codonindex : {codon : BH-adjusted pvalue}}
	psitereadcount = 0 #number of reads for which we found a valid P site
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]

	for codon in allcodons:
		totalbackgroundcounts[codon] = 0
		totalbackgroundcodondict[codon] = 0
	
	#Genecountdict was made using all "valid" reads (those that had an in-frame P site).
	#However, not all of those reads end up in codondict, because for some of them, we can't 
	#easily figure out where in the CDS the read maps (likely due to a sequencing error preventing a perfect match of the read seq somewhere within CDSseq)
	#This filtering happens within getgenestringmatch().
	#So we need to filter the reads in genecountdict so that it only contains reads present in codondict
	#Otherwise, our background codon dict would be a little off.
	print 'Filtering genecount dict...'
	readsincodondict = []
	for codonindex in codondict:
		for codon in codondict[codonindex]:
			readsincodondict += codondict[codonindex][codon]

	readsincodondict = set(readsincodondict)

	filteredgenecountdict = {} #{gene : [readIDs that map to it]}
	for gene in genecountdict:
		filteredreads = [read for read in set(genecountdict[gene]) if read in readsincodondict]
		filteredgenecountdict[gene] = filteredreads
	print 'Done!'

	#Prepare expected background dict. For every time a gene has a read, add the expected codon freq for that gene to the overall background dict
	print 'Preparing expected background codon freqs...'
	for gene in filteredgenecountdict:
		genereadcount = len(filteredgenecountdict[gene])
		psitereadcount += genereadcount
		try:
			genebgfreqs = backgroundcodonfreqs[gene]
			for codon in genebgfreqs:
				genefreq = genebgfreqs[codon]
				weightedfreq = genefreq * genereadcount
				totalbackgroundcodondict[codon] += weightedfreq
		except KeyError:
			print gene
		#Add counts of every codon to totalbackgroundcounts
		for codon in backgroundcodoncounts[gene]:
			codoncountingene = backgroundcodoncounts[gene][codon]
			weightedcodoncounts = codoncountingene * genereadcount
			totalbackgroundcounts[codon] += weightedcodoncounts

	#Normalize totalbackgroundcodondict by number of reads to get back to frequencies
	for codon in totalbackgroundcodondict:
		freq = totalbackgroundcodondict[codon] / float(psitereadcount)
		totalbackgroundcodondict[codon] = freq
	print 'Done!'

	#print totalbackgroundcodondict, sum(totalbackgroundcodondict.values())

	#prepare codoncounts, enrichments, pvalues, and adjustedpvalues dictionaries, right now they are just blank
	#They need to be of the form {codonindex : {codon : value}}
	for i in range(min(codondict.keys()), max(codondict.keys()) + 1):
		codoncounts[i] = {}
		enrichments[i] = {}
		pvalues[i] = OrderedDict()
		adjustedpvalues[i] = OrderedDict()

	#Now, for each codonindex, compare the freqs for each codon to the expected freqs in totalbackgroundcodondict
	for codonindex in codondict:
		totalxsitereads = sum(len(lst) for lst in codondict[codonindex].values()) #number of reads across all codons at this codonindex
		for codon in codondict[codonindex]:
			xsitereads = len(codondict[codonindex][codon])
			#print codon, psitereads
			othercodonxsitereads = totalxsitereads - xsitereads
			#Expected number of reads is background freq * number of total P site reads
			expectedxsitereads = totalbackgroundcodondict[codon] * totalxsitereads
			othercodonexpectedxsitereads = totalxsitereads - expectedxsitereads
			###########################################################
			####################CONTINGENCY TABLE #####################
			################RPF reads        Expected reads##############
			################_________########______________##############
			#this codon                  |
			#____________________________|_____________________##########
			#all other codons            | 
			############################################################
			contingencytable = [[xsitereads, expectedxsitereads], [othercodonxsitereads, othercodonexpectedxsitereads]]
			#Fisher exact (slow)
			#oddsratio, pvalue = stats.fisher_exact(contingencytable)
			#chi2 (faster)
			chi2, pvalue, dof, expected = stats.chi2_contingency(contingencytable)
			pvalues[codonindex][codon] = '{:.3e}'.format(pvalue)
			try:
				enrichment = (xsitereads + 1) / float(expectedxsitereads) #added pseudocount
				enrichments[codonindex][codon] = '{:.3f}'.format(log(enrichment, 2))
			except (ZeroDivisionError, ValueError):
				enrichment = 'NA'
				enrichments[codonindex][codon] = 'NA'
			codoncounts[codonindex][codon] = [xsitereads, totalxsitereads, totalbackgroundcounts[codon], sum(totalbackgroundcounts.values())]
		#print codon, psitereads, expectedpsitereads, psitereads / float(expectedpsitereads), pvalue

	#multiple hypothesis correction using BH
	#this works but we don't need it right now because we aren't using the pvalues
	uncorrectedpvals = []
	for codonindex in pvalues:
		for codon in pvalues[codonindex]:
			uncorrectedpvals.append(pvalues[codonindex][codon])

	uncorrectedpvals = [float(pval) for pval in uncorrectedpvals]
	correctedpvals = mt.multipletests(uncorrectedpvals, alpha = 0.05, method = 'fdr_bh')[1]
	#Turn into formatted strings
	correctedpvals = ['{:.3e}'.format(correctedpval) for correctedpval in correctedpvals]
	ptoadjp = dict(zip(uncorrectedpvals, correctedpvals)) #{uncorrectedpvalue : correctedpvalue}
	#populate adjustedpvalues
	for codonindex in pvalues:
		for codon in pvalues[codonindex]:
			uncorrectedpval = float(pvalues[codonindex][codon])
			correctedpval = ptoadjp[uncorrectedpval]
			adjustedpvalues[codonindex][codon] = correctedpval
	

	return enrichments, pvalues, adjustedpvalues, codoncounts


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', type = str, help = 'GFF of genome.')
    parser.add_argument('--genomefasta', type = str, help = 'Fasta of genome.')
    #parser.add_argument('--bam', type = str, help = 'Bam of RPF reads.')
    args = parser.parse_args()

    bamdir = '/beevol/home/taliaferro/data/MIT/CADFmr1KO_RFP/STARruns/bams/STARuniquemap/'
    #bamnames = ['KO-FBSrpfAAligned.sorted.bam', 'KO-FBSrpfBAligned.sorted.bam', 'KO-FBSrpfCAligned.sorted.bam', 'KO+FBSrpfAAligned.sorted.bam', 'KO+FBSrpfBAligned.sorted.bam', 'KO+FBSrpfCAligned.sorted.bam', 'WT-FBSrpfAAligned.sorted.bam', 'WT-FBSrpfBAligned.sorted.bam', 'WT-FBSrpfCAligned.sorted.bam', 'WT+FBSrpfAAligned.sorted.bam', 'WT+FBSrpfBAligned.sorted.bam', 'WT+FBSrpfCAligned.sorted.bam']
    bamnames = ['KO+FBSrpfCAligned.sorted.bam', 'WT-FBSrpfAAligned.sorted.bam', 'WT-FBSrpfBAligned.sorted.bam', 'WT-FBSrpfCAligned.sorted.bam', 'WT+FBSrpfAAligned.sorted.bam', 'WT+FBSrpfBAligned.sorted.bam', 'WT+FBSrpfCAligned.sorted.bam']
    bams = [bamdir + bam for bam in bamnames]
    print 'Loading pfd...'
    
    print 'Making CDS bed and fasta...'
    #makecdsbedandfasta(args.gff, args.genomefasta)
    print 'Getting position frame dictionary...'
    #positionframedict = getCDScoords(args.gff)
    #with open('framedict.pkl', 'wb') as outfh:
		#pickle.dump(positionframedict, outfh, protocol=pickle.HIGHEST_PROTOCOL)
	
	
    with open('framedict.pkl', 'rb') as infh:
    	positionframedict = pickle.load(infh)
    print 'Done!'

    print 'Making CDS dict...'
    longestcdsdict = makeseqdict()
    print 'Done!'


    for bam in bams:
    	print bam
	    #First filter to get the reads we will work with
    	#That's those that map to a CDS, of the right length, in frame
    	getvalidreads(bam, positionframedict, 'longestcds.bed')

    	#OK now get how many reads map to each gene, and for each read, which gene it maps to
    	print 'Getting read couts per gene...'
    	genecountdict, reads2gene = getcountspergene_tabix(os.path.basename(bam) + 'validreads.bam', 'longestcds.bed')
    	#with open('genecounts.pkl', 'wb') as outfh:
    		#pickle.dump(genecountdict, outfh, protocol=pickle.HIGHEST_PROTOCOL)
    	#with open('reads2gene.pkl', 'wb') as outfh:
    		#pickle.dump(reads2gene, outfh, protocol=pickle.HIGHEST_PROTOCOL)
    	
    	#with open('genecounts.pkl', 'rb') as infh:
    		#genecountdict = pickle.load(infh)
    	#with open('reads2gene.pkl', 'rb') as infh:
    		#reads2gene = pickle.load(infh)

    	#Now for each read, get the codons in the window surrounding the psite codon
    	print 'Iterating reads...'
    	codondict = iteratereads(os.path.basename(bam) + 'validreads.bam', positionframedict, 'longestcds.bed', reads2gene, longestcdsdict)
    	#with open('codondict.pkl', 'wb') as outfh:
    		#pickle.dump(codondict, outfh, protocol=pickle.HIGHEST_PROTOCOL)
    	#with open('codondict.pkl', 'rb') as infh:
    		#codondict = pickle.load(infh)
    	print 'Getting background codon freqs...'
    	backgroundcodonfreqs, backgroundcodoncounts = getbackgroundcodonfreqs('longestcds.fa')
    	print 'Comparing enrichments...'
    	enrichments, pvalues, adjustedpvalues, codoncounts = comparefreqs(codondict, genecountdict, backgroundcodonfreqs, backgroundcodoncounts)
    	with open('{0}.enrichments.txt'.format(os.path.basename(bam)), 'w') as outfh:
    		outfh.write(('\t').join(['codonindex', 'codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroundcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
    		for codonindex in enrichments:
    			for codon in enrichments[codonindex]:
    				#{codonindex : {codon : [number of counts in xsite, number of total counts (all codons), number of counts in background, number of total counts (all codons) in background]}}
    				sitecounts = str(codoncounts[codonindex][codon][0])
    				totalsitecounts = str(codoncounts[codonindex][codon][1])
    				backgroundcounts = str(codoncounts[codonindex][codon][2])
    				totalbackgroundcounts = str(codoncounts[codonindex][codon][3])
    				enrichment = str(enrichments[codonindex][codon])
    				pvalue = str(pvalues[codonindex][codon])
    				adjustedpvalue = str(adjustedpvalues[codonindex][codon])
    				outfh.write(('\t').join([str(codonindex), codon, sitecounts, totalsitecounts, backgroundcounts, totalbackgroundcounts, enrichment, pvalue, adjustedpvalue]) + '\n')
    
    sys.exit()
    #genecountdict, reads2gene = getcountspergene_tabix('{0}validreads.bam'.format(os.path.basename(bam)), 'longestcds.bed')

    for bam in bams:
	    print 'Calculating psites...'
	    psitedict, codondict = iteratereads(bam, positionframedict, 'longestcds.bed')
	    print 'Calculating the number of reads that map to each CDS...'
	    genecountdict = getcountspergene_tabix('{0}validreads.bam'.format(os.path.basename(bam)), 'longestcds.bed')
	    print 'Getting background codon freqs...'
	    backgroundcodonfreqs, backgroundcodoncounts = getbackgroundcodonfreqs('longestcds.fa')
	    print 'Calculating enrichments in a sites...'
	    a_enrichments, a_pvalues, a_adjustedpvalues, a_codoncounts = comparefreqs('single', asitedict, genecountdict, backgroundcodonfreqs, backgroundcodoncounts)
	    print 'Calculating enrichments in p sites...'
	    p_enrichments, p_pvalues, p_adjustedpvalues, p_codoncounts = comparefreqs('single', psitedict, genecountdict, backgroundcodonfreqs, backgroundcodoncounts)
	    print 'Calculating enrichments in e sites...'
	    e_enrichments, e_pvalues, e_adjustedpvalues, e_codoncounts = comparefreqs('single', esitedict, genecountdict, backgroundcodonfreqs, backgroundcodoncounts)
	    print 'Calculating enrichments for dicodons...'
	    di_enrichments, di_pvalues, di_adjustedpvalues, di_codoncounts = comparefreqs('double', dicodondict, genecountdict, backgrounddicodonfreqs, backgrounddicodoncounts)
	    #print 'Calculating enrichments for tricodons...'
	    #tri_enrichments, tri_pvalues, tri_adjustedpvalues, tri_codoncounts = comparefreqs(tricodondict, genecountdict, backgroundtricodonfreqs, backgroundtricodoncounts)

	    #Write output
	    with open('{0}.asiteenrichments.txt'.format(os.path.basename(bam)), 'w') as outfh:
	    	outfh.write(('\t').join(['codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroundcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
	    	for codon in a_enrichments:
	    		outfh.write(('\t').join([codon, str(a_codoncounts[codon][0]), str(a_codoncounts[codon][1]), str(a_codoncounts[codon][2]), str(a_codoncounts[codon][3]), a_enrichments[codon], a_pvalues[codon], a_adjustedpvalues[codon]]) + '\n')

	    with open('{0}.psiteenrichments.txt'.format(os.path.basename(bam)), 'w') as outfh:
	    	outfh.write(('\t').join(['codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroundcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
	    	for codon in p_enrichments:
	    		outfh.write(('\t').join([codon, str(p_codoncounts[codon][0]), str(p_codoncounts[codon][1]), str(p_codoncounts[codon][2]), str(p_codoncounts[codon][3]), p_enrichments[codon], p_pvalues[codon], p_adjustedpvalues[codon]]) + '\n')

	    with open('{0}.esiteenrichments.txt'.format(os.path.basename(bam)), 'w') as outfh:
	    	outfh.write(('\t').join(['codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroundcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
	    	for codon in e_enrichments:
	    		outfh.write(('\t').join([codon, str(e_codoncounts[codon][0]), str(e_codoncounts[codon][1]), str(e_codoncounts[codon][2]), str(e_codoncounts[codon][3]), e_enrichments[codon], e_pvalues[codon], e_adjustedpvalues[codon]]) + '\n')

	    with open('{0}.disiteenrichments.txt'.format(os.path.basename(bam)), 'w') as outfh:
	    	outfh.write(('\t').join(['codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroundcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
	    	for codon in di_enrichments:
	    		outfh.write(('\t').join([codon, str(di_codoncounts[codon][0]), str(di_codoncounts[codon][1]), str(di_codoncounts[codon][2]), str(di_codoncounts[codon][3]), di_enrichments[codon], di_pvalues[codon], di_adjustedpvalues[codon]]) + '\n')

	    #with open('{0}.trisiteenrichments.txt'.format(os.path.basename(args.bam)), 'w') as outfh:
	    	#outfh.write(('\t').join(['codon', 'sitecounts', 'totalsitecounts', 'backgroundcounts', 'totalbackgroudcounts', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
	    	#for codon in tri_enrichments:
	    		#outfh.write(('\t').join([codon, str(tri_codoncounts[codon][0]), str(tri_codoncounts[codon][1]), str(tri_codoncounts[codon][2]), str(tri_codoncounts[codon][3]), str(tri_enrichments[codon]), str(tri_pvalues[codon]), str(tri_adjustedpvalues[codon])]) + '\n')

