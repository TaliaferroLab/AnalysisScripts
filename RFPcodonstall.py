import gffutils
import os
import sys
import pysam
import subprocess
import argparse
import pickle
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

	with open('framedict.pkl', 'wb') as outfh:
		pickle.dump(framedict, outfh)

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
			if codon not in allcodons:
				continue
			backgroundcodoncounts[record.id][codon] += 1

	#Turn counts into frequencies
	for gene in backgroundcodoncounts:
		backgroundcodonfreqs[gene] = {}
		totalcodoncounts = sum(backgroundcodoncounts[gene].values())
		for codon in backgroundcodoncounts[gene]:
			freq = backgroundcodoncounts[gene][codon] / float(totalcodoncounts)
			backgroundcodonfreqs[gene][codon] = freq


	return backgroundcodonfreqs


#Given a read, get its P site codon.
def getcodon(read, positionframedict):
	if read.is_reverse:
		strand = '-'
	elif not read.is_reverse:
		strand = '+'
	chrm = read.reference_name
	#No mitochondrial ribosome stuff
	if chrm == 'chrM':
		return None, None, None
	if strand == '+':
		readstart = read.reference_start + 1
	elif strand == '-':
		readstart = read.reference_end

	#Offset readstart by 13 to get P site
	refpositions = read.get_reference_positions(full_length = True)
	#Using full_length = True will allow soft-clipped nts to have placeholders of None in this list
	#We still want them as placeholders so that the offset will be correct
	refpositions = [pos + 1 if type(pos) == int else pos for pos in refpositions] #make 1 based
	if read.query_length == 28 or read.query_length == 29 or read.query_length == 30:
	#if read.query_length == 20 or read.query_length == 21:
		if strand == '+':
			refpos = refpositions[13]
		elif strand == '-':
			refpositions.reverse()
			refpos = refpositions[13]

	#Not all readstarts are going to be in positionframedict.
	#Not all reads map to CDS.
	try:
		frame = positionframedict[chrm][strand][refpos]
		#If the P site for this read is in frame 0, get the codon in the P site
		if frame == 0:
			if strand == '+':
				readseq = read.query_sequence
			elif strand == '-':
				readseq = read.get_forward_sequence()
			esitecodon = readseq[10:13]
			psitecodon = readseq[13:16]
			asitecodon = readseq[16:19]
			#print chrm, readstart, strand, read.query_sequence, codon
		elif frame != 0:
			psitecodon, asitecodon, esitecodon = None, None, None
	except KeyError:
		psitecodon, asitecodon, esitecodon = None, None, None

	return psitecodon, asitecodon, esitecodon

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

#Given a set of gene coordinates (or CDS coords) and a bam, get the reads that map to each gene
def getcountspergene(bam, cdsbed):
	genecounts = {} #{gene : [readIDs that map to it]}
	bambedtool = pybedtools.BedTool(bam)
	with open(cdsbed, 'r') as cdsbedfh:
		#Go through each line of the bed (gene) individually, making a new one-line bed file
		#Then intersect that with the bam
		genecounter = 0
		for line in cdsbedfh:
			genecounter +=1
			if genecounter % 100 == 0:
				print 'Gene {0}...'.format(genecounter)
			line = line.strip().split('\t')
			genename = line[3]
			genebed = (' ').join(line)

			genecounts[genename] = []
			genebedtool = pybedtools.BedTool(genebed, from_string = True)
			readmatches = bambedtool.intersect(genebedtool, s = True, bed = True)
			for readmatch in readmatches:
				readid = str(readmatch).split('\t')[3]
				genecounts[genename].append(readid)

	with open('{0}genecounts.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(genecounts, outfh)

	return genecounts

def getcountspergene_tabix(bam, cdsbed):
	genecounts = {} #{gene : [readIDs that map to it]}

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
				genecounts[genename].append(readid)

	with open('{0}genecounts.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(genecounts, outfh)

	os.remove('temp.sam.gz')
	os.remove('temp.sam.gz.tbi')

	return genecounts



#Given the frame of all CDS positions and the bam, get the starting position of all reads
#and then get their frames
#For each gene, record the number of times each codon is in the P site
def iteratereads(bam, positionframedict, cdsbed):
	Psitedict = {} #{codon in P site : list of readIDs that have this as a psite codon}
	Asitedict = {} #{codon in A site : list of readIDs that have this as a asite codon}
	Esitedict = {} #{codon in E site : list of readIDs that have this as a esite codon}
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]
	for codon in allcodons:
		Psitedict[codon] = []
		Asitedict[codon] = []
		Esitedict[codon] = []


	cdsbedtool = pybedtools.BedTool(cdsbed)
	readcounter = 0
	with pysam.AlignmentFile(bam, 'rb') as infh, pysam.AlignmentFile('{0}readswithPsitecodons.bam'.format(os.path.basename(bam)), 'wb', template = pysam.AlignmentFile(bam, 'rb')) as outfh:
		for read in infh.fetch(until_eof = True):
			readcounter +=1
			if readcounter % 1000000 == 0:
				print 'Read {0}...'.format(readcounter)
			#Only consider reads of lengths 30-35
			#if read.query_length not in [30, 31, 32, 33, 34, 35]:
			#if read.query_length not in [28, 29, 30]:
			if read.query_length not in [20, 21]:
				continue

			#For non-RFP reads, the first read of mate pair is on the wrong strand
			#Also, only consider unique mappers
			if read.is_read1 or read.get_tag('NH') > 1 or read.is_secondary:
				continue

			psitecodon, asitecodon, esitecodon = getcodon(read, positionframedict)
			#gene = getgene(read, cdsbedtool)

			#If this read maps to one longestCDS and is in the 0 frame (meaning it returns something from getcodon), we will deal with it
			#That means neither codon nor gene are 'None'
			#if codon and gene:
			if psitecodon in allcodons and asitecodon in allcodons and esitecodon in allcodons:
				#print codon, gene
				#if gene not in Psitedict:
					#Psitedict[gene] = {}
					#for codon in allcodons:
						#Psitedict[gene][codon] = 0

				#Psitedict[gene][codon] += 1
				Psitedict[psitecodon].append(read.query_name)
				Asitedict[asitecodon].append(read.query_name)
				Esitedict[esitecodon].append(read.query_name)
				outfh.write(read)

	with open('{0}psite.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		print 'YES'
		pickle.dump(Psitedict, outfh)
	with open('{0}asite.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(Asitedict, outfh)
	with open('{0}esite.pkl'.format(os.path.basename(bam)), 'wb') as outfh:
		pickle.dump(Esitedict, outfh)
	return Psitedict, Asitedict, Esitedict

#Compare the frequencies for a codon being in the P site to expected background across all genes
#These frequencies are weighted by the number of times a gene had a Psite read in it
#Each gene has its own bg freq (that's in backgroundcodonfreqs).
#These freqs are multiplied by the number of Psite reads in that gene (then at the end divided by the total number of P site reads to get back to freqs).
def comparefreqs(xsitedict, genecountdict, backgroundcodonfreqs):
	#Psitedict = {} #{codon in P site : list of readIDs that have this as a psite codon}
	#genecountdict = {} #{gene : [readIDs that map to it]}
	#backgroundcodonfreqs = {} #{geneid : {codon : freq in that gene}}
	totalbackgroundcodondict = {} #{codon: final expected bg freq (all of these should add up to 1)}
	enrichments = {} #{codon : log2 enrichment over background}
	pvalues = OrderedDict() #{codon : pvalue}
	adjustedpvalues = OrderedDict() #{codon : BH-adjusted pvalue}
	psitereadcount = 0 #number of reads for which we found a valid P site
	allcodons = [''.join(codon) for codon in product(['A', 'C', 'T', 'G'], repeat = 3)]

	for codon in allcodons:
			totalbackgroundcodondict[codon] = 0
	
	#Prepare expected background dict. For every time a gene has a read, add the expected codon freq for that gene to the overall background dict
	for gene in genecountdict:
		genereadcount = len(genecountdict[gene])
		psitereadcount += genereadcount
		try:
			genebgfreqs = backgroundcodonfreqs[gene]
			for codon in genebgfreqs:
				genefreq = genebgfreqs[codon]
				weightedfreq = genefreq * genereadcount
				totalbackgroundcodondict[codon] += weightedfreq
		except KeyError:
			print gene

	#Normalize totalbackgroundcodondict by number of reads to get back to frequencies
	for codon in totalbackgroundcodondict:
		freq = totalbackgroundcodondict[codon] / float(psitereadcount)
		totalbackgroundcodondict[codon] = freq

	#print totalbackgroundcodondict, sum(totalbackgroundcodondict.values())

	#Now compare the freqs in Psitedict to the expected freqs in totalbackgroundcodondict
	totalxsitereads = sum(len(lst) for lst in xsitedict.values())
	for codon in xsitedict:
		xsitereads = len(xsitedict[codon])
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
		oddsratio, pvalue = stats.fisher_exact(contingencytable)
		pvalues[codon] = pvalue
		enrichment = xsitereads / float(expectedxsitereads)
		try:
			enrichments[codon] = log(enrichment, 2)
		except ValueError: #log of 0
			enrichments[codon] = 'NA'
		#print codon, psitereads, expectedpsitereads, psitereads / float(expectedpsitereads), pvalue

	#multiple hypothesis correction using BH
	#this works but we don't need it right now because we aren't using the pvalues
	
	uncorrectedpvals = pvalues.values()
	correctedpvals = mt.multipletests(uncorrectedpvals, alpha = 0.05, method = 'fdr_bh')[1]
	adjustedpvalues = dict(zip(pvalues.keys(), correctedpvals))
	for codon in pvalues:
		#print codon, pvalues[codon], adjustedpvalues[codon], enrichments[codon]
		pass
	

	return enrichments, pvalues, adjustedpvalues

#Combine sitedicts and genecountdicts across samples
def combinedicts():
	asites = {} #{codon : [list of reads with asite codon]}
	psites = {} #{codon : [list of reads with psite codon]}
	genecounts = {} #{gene : [list of reads in that gene]}

	for sample in ['E301', 'E302', 'E601', 'E602', 'E901', 'E902', 'P301', 'P302', 'P601', 'P602', 'P901', 'P902']:
		print sample
		asitedict = '{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bamasite.pkl'.format(sample)
		psitedict = '{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bampsite.pkl'.format(sample)
		genecountdict = '{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bamreadswithPsitecodons.bamgenecounts.pkl'.format(sample)

		with open(asitedict, 'rb') as asitefh, open(psitedict, 'rb') as psitefh, open(genecountdict, 'rb') as genecountfh:
			asitedict = pickle.load(asitefh)
			psitedict = pickle.load(psitefh)
			genecountdict = pickle.load(genecountfh)

			for codon in asitedict:
				if codon not in asites:
					asites[codon] = []
				asites[codon] += asitedict[codon]

			for codon in psitedict:
				if codon not in psites:
					psites[codon] = []
				psites[codon] += psitedict[codon]

			for gene in genecountdict:
				if gene not in genecounts:
					genecounts[gene] = []
				genecounts[gene] += genecountdict[gene]

	#Remove empty (nonexpressed) genes
	genecounts_exp = {k: v for k, v in genecounts.items() if v}

	with open('allasites.pkl', 'wb') as asitefh, open('allpsites.pkl', 'wb') as psitefh, open('allgenecounts.pkl', 'wb') as genecountfh:
		pickle.dump(asites, asitefh)
		pickle.dump(psites, psitefh)
		pickle.dump(genecounts_exp, genecountfh)

#Suppose we only wanted to look at these enrichments across a subset of genes and not all of them.
#We can first look in the genecountdict to see which reads map to those genes
#Then subset the psitedict and asitedict so that only those reads are in them.
def subsetgenes(genelist, genecountdict, psitedict, asitedict):
	#Subset genecountdict based on genelist
	subgenecountdict = {gene: genecountdict[gene] for gene in genelist}
	#Get all the reads that map to these genes
	subreads = subgenecountdict.values() #this is a list of lists so we need to flatten it
	subreads = [read for gene in subreads for read in gene]
	subreads = set(subreads) #for easy intersection

	#Now go through each codon intersecting each codons reads with subreads
	subpsitedict = {} #{codon : subsetted reads that have that codon as psite}
	subasitedict = {} #{codon : subsetted reads that have that codon as asite}
	for codon in psitedict:
		codonreads = set(psitedict[codon])
		subcodonreads = codonreads.intersection(subreads)
		subpsitedict[codon] = list(subcodonreads)
	
	for codon in asitedict:
		codonreads = set(asitedict[codon])
		subcodonreads = codonreads.intersection(subreads)
		subasitedict[codon] = list(subcodonreads)

	return subgenecountdict, subpsitedict, subasitedict

def getsubsetenrichments(genecountdict, asitedict, psitedict, backgroundcodonfreqs):
	allgenes = genecountdict.keys()
	with open('subsetenrichments.txt', 'w') as outfh:
		outfh.write(('\t').join(['sample', 'codon', 'site', 'enrichment']) + '\n')
		for i in range(500):
			print i + 1
			samplegenes = sample(allgenes, 100)
			subgenecountdict, subpsitedict, subasitedict = subsetgenes(samplegenes, genecountdict, psitedict, asitedict)
			asiteenrichments = comparefreqs(subasitedict, subgenecountdict, backgroundcodonfreqs)
			psiteenrichments = comparefreqs(subpsitedict, subgenecountdict, backgroundcodonfreqs)
			for codon in asiteenrichments:
				outfh.write(('\t').join([str(i + 1), codon, 'a', str(asiteenrichments[codon])]) + '\n')
			for codon in psiteenrichments:
				outfh.write(('\t').join([str(i + 1), codon, 'p', str(psiteenrichments[codon])]) + '\n')





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', type = str, help = 'GFF of genome.')
    parser.add_argument('--genomefasta', type = str, help = 'Fasta of genome.')
    parser.add_argument('--bam', type = str, help = 'Bam of RPF reads.')
    args = parser.parse_args()

    print 'Making CDS bed and fasta...'
    makecdsbedandfasta(args.gff, args.genomefasta)
    print 'Getting position frame dictionary...'
    positionframedict = getCDScoords(args.gff)
    print 'Calculating psites and asites...'
    psitedict, asitedict, esitedict = iteratereads(args.bam, positionframedict, 'longestcds.bed')
    print 'Calculating the number of reads that map to each CDS...'
    genecountdict = getcountspergene_tabix('{0}readswithPsitecodons.bam'.format(os.path.basename(args.bam)), 'longestcds.bed')
    print 'Getting background codon freqs...'
    backgroundcodonfreqs = getbackgroundcodonfreqs('longestcds.fa')
    print 'Calculating enrichments in a sites...'
    a_enrichments, a_pvalues, a_adjustedpvalues = comparefreqs(asitedict, genecountdict, backgroundcodonfreqs)
    print 'Calculating enrichments in p sites...'
    p_enrichments, p_pvalues, p_adjustedpvalues = comparefreqs(psitedict, genecountdict, backgroundcodonfreqs)
    print 'Calculating enrichments in e sites...'
    e_enrichments, e_pvalues, e_adjustedpvalues = comparefreqs(esitedict, genecountdict, backgroundcodonfreqs)

    #Write output
    with open('{0}.asiteenrichments.txt'.format(os.path.basename(args.bam)), 'w') as outfh:
    	outfh.write(('\t').join(['codon', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
    	for codon in a_enrichments:
    		outfh.write(('\t').join([codon, str(a_enrichments[codon]), str(a_pvalues[codon]), str(a_adjustedpvalues[codon])]) + '\n')

    with open('{0}.psiteenrichments.txt'.format(os.path.basename(args.bam)), 'w') as outfh:
    	outfh.write(('\t').join(['codon', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
    	for codon in p_enrichments:
    		outfh.write(('\t').join([codon, str(p_enrichments[codon]), str(p_pvalues[codon]), str(p_adjustedpvalues[codon])]) + '\n')

    with open('{0}.esiteenrichments.txt'.format(os.path.basename(args.bam)), 'w') as outfh:
    	outfh.write(('\t').join(['codon', 'enrichment', 'pvalue', 'adjustedpvalue']) + '\n')
    	for codon in e_enrichments:
    		outfh.write(('\t').join([codon, str(e_enrichments[codon]), str(e_pvalues[codon]), str(e_adjustedpvalues[codon])]) + '\n')

    '''
    for sample in ['E902', 'P901', 'P302', 'E302', 'P301', 'P602', 'P601']:
    	with open('{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bampsite.pkl'.format(sample), 'rb') as infh:
    		psitedict = pickle.load(infh)
   		with open('{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bamasite.pkl'.format(sample), 'rb') as infh:
   			asitedict = pickle.load(infh)
    	with open('{0}.UniqueMap.RFP.Aligned.sortedByCoord.out.bamreadswithPsitecodons.bamgenecounts.pkl'.format(sample), 'rb') as infh:
    		genecountdict = pickle.load(infh)

    genelist = ['ENSG00000162772', 'ENSG00000134107', 'ENSG00000118523', 'ENSG00000142871', 'ENSG00000120129', 'ENSG00000120738',
    'ENSG00000063046', 'ENSG00000160888', 'ENSG00000177606', 'ENSG00000171223', 'ENSG00000155090', 'ENSG00000123358', 'ENSG00000143878', 
    'ENSG00000142676', 'ENSG00000148303', 'ENSG00000141232', 'ENSG00000185650']
    subgenecountdict, subpsitedict, subasitedict = subsetgenes(genelist, genecountdict, psitedict, asitedict)

    	penrichments = comparefreqs(psitedict, genecountdict, backgroundcodonfreqs)
    	aenrichments = comparefreqs(asitedict, genecountdict, backgroundcodonfreqs)
    	with open('{0}penrichments.txt'.format(sample), 'w') as outfh:
    		outfh.write(('\t').join(['codon', 'enrichment']) + '\n')
    		for codon in penrichments:
    			outfh.write(('\t').join([codon, str(penrichments[codon])]) + '\n')

    	with open('{0}aenrichments.txt'.format(sample), 'w') as outfh:
    		outfh.write(('\t').join(['codon', 'enrichment']) + '\n')
    		for codon in aenrichments:
    			outfh.write(('\t').join([codon, str(aenrichments[codon])]) + '\n')

    #Combine asitedicts and psitedicts and genecountdicts
    #combinedicts()
    #Run subsets
    with open('allgenecounts.pkl', 'rb') as infh:
    	print 'Loading gene counts...'
    	allgenecounts = pickle.load(infh)
    with open('allasites.pkl', 'rb') as infh:
    	print 'Loading a sites...'
    	allasites = pickle.load(infh)
    with open('allpsites.pkl', 'rb') as infh:
    	print 'Loading p sites...'
    	allpsites = pickle.load(infh)

    getsubsetenrichments(allgenecounts, allasites, allpsites, backgroundcodonfreqs)
	'''
