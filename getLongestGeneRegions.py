#Given a genome annotation (ensembl preferably) and a genome sequence, get the longest
#3' UTR, 5' UTR, and CDS coords for every gene in the annotation.

import gffutils
import os
from Bio import SeqIO
import argparse
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import gzip
import sys

#########################
#########################
#Get longest CDS of gene#
#########################
#########################
def getCDScoords(gff, ens2short, outputgff):
	allCDScoords = {} #{ENSGENE_chrm_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	geneswithcodingtranscript = 0
	e2sdict = {} #{ENSGene : shortname}

	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[2]
	infh.close()

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

	print 'Looked through {0} genes. {1} of them had non-NMD coding transcripts. Found longest CDS sequences for {2} of them.'.format(genecount, geneswithcodingtranscript, len(allCDScoords))
	with open(outputgff, 'w') as f:
		geneswithshortnames = 0
		for gene in allCDScoords:
			exoncounter = 0
			ensid = gene.split('_')[0].replace('gene:', '') #in allCDScoords genenames are like 'gene:ENSMUSG00000038375_chr2_+'
			if ensid not in e2sdict:
				shortname = ensid
			else:
				shortname = e2sdict[ensid]

			chrm = gene.split('_')[1]
			strand = gene.split('_')[2]
			genestart, geneend = geneboundaries[ensid][0], geneboundaries[ensid][1]
			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(ensid, shortname, ensid)
			f.write(('\t').join([chrm, 'longestCDS', 'gene', str(genestart), str(geneend), '.', strand, '.', IDline]) + '\n')
			CDSstart = allCDScoords[gene][0][0]
			CDSstop = allCDScoords[gene][-1][1]
			IDline = 'ID=CDS:{0};Parent=gene:{1};gene_id={2}'.format(ensid, ensid, shortname)
			f.write(('\t').join([chrm, 'longestCDS', 'CDS', str(CDSstart), str(CDSstop), '.', strand, '.', IDline]) + '\n')
			for CDSexon in allCDScoords[gene]:
				exoncounter +=1
				exonstart = CDSexon[0]
				exonstop = CDSexon[1]
				IDline = 'ID=exon:{0}.cdsexon{1};Parent=CDS:{2}'.format(ensid, exoncounter, ensid)
				f.write(('\t').join([chrm, 'longestCDS', 'exon', str(exonstart), str(exonstop), '.', strand, '.', IDline]) + '\n')




	return allCDScoords


########################
#Turn longest CDS coords into sequence#
########################

def getCDSSequences(allCDScoords, genomefasta, ens2short, outputfasta):
	#ens2short for mm10 is Ensembl_to_genename.txt
	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	CDSseqs = {} #{genename : CDSseq}
	e2sdict = {} #{ENSGene : shortname}
	chrmswithoutseq = [] #Chromsome names that are in allCDScoords but that don't have a fasta entry in genomefasta

	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[2]
			#e2sdict[line[1]] = line[0]
	infh.close()

	for gene in allCDScoords:
		CDSseq = ''
		genename = gene.split('_')[0].replace('gene:','') #in allCDScoords genenames are like 'gene:ENSMUSG00000038375_chr2_+'
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]

		#Is this chromosome in genomefasta?
		if chrm not in seq_dict: 
			if chrm not in chrmswithoutseq:
				print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
				chrmswithoutseq.append(chrm)
			continue

		for coords in allCDScoords[gene]:
			start = coords[0]
			end = coords[1]
			if strand == '+':
				exonseq = seq_dict[chrm].seq[start-1:end].upper()
				CDSseq += exonseq
			elif strand == '-':
				exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
				newseq = exonseq + CDSseq
				CDSseq = newseq

		if genename in e2sdict:
			shortname = e2sdict[genename]
			CDSseqs[genename + '_' + shortname] = str(CDSseq)

		elif genename not in e2sdict:
			shortname = genename
			CDSseqs[genename + '_' + shortname] = str(CDSseq)

	with open(outputfasta, 'w') as f:
		for CDS in CDSseqs:
			f.write('>' + CDS + '\n' + str(CDSseqs[CDS]) + '\n')

	return CDSseqs


######################
#Get longest 3pUTR for every gene#
######################

def getlongest3putrcoords(gff, ens2short, outputgff):
	all3pUTRcoords = defaultdict(list) # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	longest3pUTRcoords = {} # {ENSGENE_chrm_strand : [[exon1start, exon1end], [exon2start, exon2end]]}
	e2sdict = {} #{ENSGene : shortname}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	no3pUTRcount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = gff + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Get ensembl ID / gene short name relationships
	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			#e2sdict[line[0]] = line[2]
			e2sdict[line[0]] = line[2]
	infh.close()

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
				all3pUTRcoords[genename].append(UTRexoncoords)

	#Now get longest 3' UTR for each gene
	for gene in all3pUTRcoords:
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
		longest3pUTRcoords[gene] = all3pUTRcoords[gene][indexofmax]

		#Turn single exon utrs that are [exonstart, exonstop] into [[exonstart, exonstop]]
		#so that they can be parsed like multiexon utrs
		if len(all3pUTRcoords[gene][indexofmax]) == 2 and type(all3pUTRcoords[gene][indexofmax][0]) == int:
			longest3pUTRcoords[gene] = [all3pUTRcoords[gene][indexofmax]]

	with open(outputgff, 'w') as f:
		#longest3pUTRcoords = {} # {ENSGENE_chrm_strand : [[exon1start, exon1end], [exon2start, exon2end]]}
		for gene in longest3pUTRcoords:
			ensid = gene.split('_')[0].replace('gene:', '')
			chrm = gene.split('_')[1]
			strand = gene.split('_')[2]
			if ensid in e2sdict:
				shortname = e2sdict[ensid]
			else:
				print 'WARNING: no ENSID found for {0}.'.format(ensid)
				shortname = ensid
			#shortname = gene.split('_')[0].split(':')[1]
			#chrm = gene.split('_')[1]
			#strand = gene.split('_')[2]
			#if shortname in e2sdict:
				#ensid = e2sdict[shortname]
			#else:
				#print 'WARNING: no ENSID found for {0}.'.format(shortname)
				#ensid = shortname
			genestart = geneboundaries[ensid][0]
			geneend = geneboundaries[ensid][1]
			#genestart = geneboundaries[shortname][0]
			#geneend = geneboundaries[shortname][1]
			geneIDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(ensid, shortname, ensid)
			UTRIDline = 'ID=UTR3:{0};Parent=gene:{1};gene_id={2}'.format(ensid, ensid, shortname)
			
			if type(longest3pUTRcoords[gene][0]) == list:
				UTRstart = longest3pUTRcoords[gene][0][0]
				UTRend = longest3pUTRcoords[gene][-1][1]
				f.write(('\t').join([chrm, 'longest3UTR', 'gene', str(genestart), str(geneend), '.', strand, '.', geneIDline]) + '\n')
				f.write(('\t').join([chrm, 'longest3UTR', 'UTR3', str(UTRstart), str(UTRend), '.', strand, '.', UTRIDline]) + '\n')
				exoncounter = 0
				for exon in longest3pUTRcoords[gene]:
					exoncounter +=1
					exonstart = exon[0]
					exonend = exon[1]
					exonIDline = 'ID=exon:{0}.utr3exon{1};Parent=UTR3:{2}'.format(ensid, exoncounter, ensid)
					f.write(('\t').join([chrm, 'longest3UTR', 'exon', str(exonstart), str(exonend), '.', strand, '.', exonIDline]) + '\n')

			else:
				print 'ERROR!'

	return longest3pUTRcoords


def getlongest5putrcoords(gff, ens2short, outputgff):
	all5pUTRcoords = defaultdict(list) # {ENSGENE_chrm_strand : [[UTR1start, UTR1end], [UTR2start, UTR2end], [[UTR3exon1start, UTR3exon1end], [UTR3exon2start, UTR3exon2end]]]}
	longest5pUTRcoords = {} # {ENSGENE_chrm_strand : [[exon1start, exon1end], [exon2start, exon2end]]}
	e2sdict = {} #{ENSGene : shortname}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	no3pUTRcount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = gff + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Get ensembl ID / gene short name relationships
	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[2]
	infh.close()

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
			for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
				exoncoords.append([exon.start, exon.end])
			for CDSexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
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
					all5pUTRcoords[genename].append([UTRcoords[0], UTRcoords[1]])

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
				#if strand == '-':
					#UTRexoncoords.reverse()
				all5pUTRcoords[genename].append(UTRexoncoords)

	#Now get longest 5' UTR for each gene
	for gene in all5pUTRcoords:
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
		longest5pUTRcoords[gene] = all5pUTRcoords[gene][indexofmax]
		if len(all5pUTRcoords[gene][indexofmax]) == 2 and type(all5pUTRcoords[gene][indexofmax][0]) == int:
			longest5pUTRcoords[gene] = [all5pUTRcoords[gene][indexofmax]]

	with open(outputgff, 'w') as f:
		#longest5pUTRcoords = {} # {ENSGENE_chrm_strand : [[exon1start, exon1end], [exon2start, exon2end]]}
		for gene in longest5pUTRcoords:
			ensid = gene.split('_')[0].replace('gene:', '')
			chrm = gene.split('_')[1]
			strand = gene.split('_')[2]
			if ensid not in e2sdict:
				shortname = ensid
			else:
				shortname = e2sdict[ensid]
			genestart = geneboundaries[ensid][0]
			geneend = geneboundaries[ensid][1]
			geneIDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(ensid, shortname, ensid)
			UTRIDline = 'ID=UTR5:{0};Parent=gene:{1};gene_id={2}'.format(ensid, ensid, shortname)

			if type(longest5pUTRcoords[gene][0]) == list:
				UTRstart = longest5pUTRcoords[gene][0][0]
				UTRend = longest5pUTRcoords[gene][-1][1]
				f.write(('\t').join([chrm, 'longest5UTR', 'gene', str(genestart), str(geneend), '.', strand, '.', geneIDline]) + '\n')
				f.write(('\t').join([chrm, 'longest5UTR', 'UTR5', str(UTRstart), str(UTRend), '.', strand, '.', UTRIDline]) + '\n')
				exoncounter = 0
				for exon in longest5pUTRcoords[gene]:
					exoncounter +=1
					exonstart = exon[0]
					exonend = exon[1]
					exonIDline = 'ID=exon:{0}.utr5exon{1};Parent=UTR5:{2}'.format(ensid, exoncounter, ensid)
					f.write(('\t').join([chrm, 'longest5UTR', 'exon', str(exonstart), str(exonend), '.', strand, '.', exonIDline]) + '\n')

			else:
				print 'ERROR!'

	return longest5pUTRcoords





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mm10gff', type = str, help = 'mm10 annotation in gff format.  Usually Mus_musculus.GRCm38.83.gff3.')
	parser.add_argument('--mm10genome', type = str, help = 'mm10 genome in fasta format.  Usually mm10.fasta.gz')
	parser.add_argument('--ens2short', type = str, help = 'File containing ensembl ID to gene short name relations. Usually Ensembl_to_genename.txt')
	parser.add_argument('--region', type = str, choices = ['CDS', 'UTR5', 'UTR3'], help = 'Gene region to go after.')
	parser.add_argument('--outputgff', type = str, help = 'Gene regions output in gff format.')
	parser.add_argument('--outputfasta', type = str, help = 'Sequence of regions in fasta format.')
	args = parser.parse_args()

	if args.region == 'CDS':
		allCDScoords = getCDScoords(args.mm10gff, args.ens2short, args.outputgff)
		getCDSSequences(allCDScoords, args.mm10genome, args.ens2short, args.outputfasta)

	elif args.region == 'UTR3':
		longest3putrcoords = getlongest3putrcoords(args.mm10gff, args.ens2short, args.outputgff)
		getCDSSequences(longest3putrcoords, args.mm10genome, args.ens2short, args.outputfasta)

	elif args.region == 'UTR5':
		longest5putrcoords = getlongest5putrcoords(args.mm10gff, args.ens2short, args.outputgff)
		getCDSSequences(longest5putrcoords, args.mm10genome, args.ens2short, args.outputfasta)







