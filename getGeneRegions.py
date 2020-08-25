#Given a genome annotation (ensembl preferably) and a genome sequence, and a list of transcript
#IDs that you are interested in, get coords and sequences for the 5' UTR, 3' UTR, and CDS
#regions of those transcripts

import gffutils
import os
from Bio import SeqIO
import argparse
import gzip
from operator import itemgetter
from itertools import groupby

def getCDScoords(gff, ens2short, txs, outputgff):
	txCDScoords = {} #{ENSMUST_chrm_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	tx2gene = {} # {ENSMUST : ENSMUSG}
	txchrmstrand = {} #{ENSMUST: [chrm, strand]}
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
		geneID = str(gene.id).replace('gene:', '')
		geneboundaries[geneID] = [gene.start, gene.end]
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start'):
			#If this transcript has no coding exons, skip it:
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				continue

			transcriptID = str(transcript.id).replace('transcript:', '').split('.')[0]
			if transcriptID in txs:
				tx2gene[transcriptID] = geneID
				#txCDScoords[transcriptID + '_' + chrm + '_' + strand] = []
				txCDScoords[transcriptID] = []
				txchrmstrand[transcriptID] = [transcript.chrom, transcript.strand]
				for codingexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
					#txCDScoords[transcriptID + '_' + chrm + '_' + strand].append([codingexon.start, codingexon.end])
					txCDScoords[transcriptID].append([codingexon.start, codingexon.end])

	print 'Looked through {0} genes for {1} transcripts. Found coding regions for {2} of them.'.format(genecount, len(txs), len(txCDScoords))

	with open(outputgff, 'w') as f:
		for transcript in txCDScoords:
			exoncounter = 0
			transcriptID = transcript.split('_')[0]
			chrm = txchrmstrand[transcriptID][0]
			strand = txchrmstrand[transcriptID][1]
			geneID = tx2gene[transcriptID]
			geneshortname = e2sdict[geneID.split('.')[0]]
			genestart, geneend = geneboundaries[geneID][0], geneboundaries[geneID][1]
			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(geneID, geneshortname, geneID)
			f.write(('\t').join([chrm, 'gene', 'gene', str(genestart), str(geneend), '.', strand, '.', IDline]) + '\n')
			CDSstart = txCDScoords[transcript][0][0]
			CDSstop = txCDScoords[transcript][-1][1]
			IDline = 'ID=CDS:{0};Parent={1};gene_id={2}'.format(transcriptID, geneID, geneID)
			f.write(('\t').join([chrm, 'CDS', 'CDS', str(CDSstart), str(CDSstop), '.', strand, '.', IDline]) + '\n')
			for CDSexon in txCDScoords[transcript]:
				exoncounter +=1
				exonstart = CDSexon[0]
				exonstop = CDSexon[1]
				IDline = 'ID=exon:{0}.cdsexon{1};Parent=CDS:{2}'.format(transcriptID, exoncounter, transcriptID)
				f.write(('\t').join([chrm, 'exon', 'exon', str(exonstart), str(exonstop), '.', strand, '.', IDline]) + '\n')

	return txCDScoords, tx2gene, txchrmstrand


def get3UTRcoords(gff, ens2short, txs, outputgff):
	txUTRcoords = {} #{ENSMUST_chrm_strand : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
	txchrmstrand = {} #{ENSMUST: [chrm, strand]}
	genechrmstrand = {} #{ENSMUSG : [chrm , strand]}
	gene2tx = {} # {ENSMUSG : [ENSMUST]}
	tx2gene = {} # {ENSMUST : ENSMSUG}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	e2sdict = {} #{ENSGene : shortname}

	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[2]
	infh.close()

	a = []
	for tx in txs:
		a.append(tx.split('.')[0])
	txs = a

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
		geneID = str(gene.id).replace('gene:', '')
		geneboundaries[geneID] = [gene.start, gene.end]
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start'):
			#If this transcript has no coding exons, skip it:
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				continue

			transcriptID = str(transcript.id).replace('transcript:', '').split('.')[0]
			
			if transcriptID in txs:

				if geneID not in genechrmstrand:
					genechrmstrand[geneID] = [gene.chrom, gene.strand]

				if transcriptID not in txchrmstrand:
					txchrmstrand[transcriptID] = [transcript.chrom, transcript.strand]

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
					#If the transcript ends right where the CDS ends, then there's no UTR
					if CDSend == transcript.end:
						continue
					UTR3start = CDSend + 1
					UTRcoords = [UTR3start, transcript.end]
				elif transcript.strand == '-':
					CDSend = min(CDScoords, key = itemgetter(0))[0]
					#If the transcript ends right where the CDS ends, then there's no UTR
					if CDSend == transcript.start:
						continue
					UTR3start = CDSend - 1
					UTRcoords = [transcript.start, UTR3start]

				#Check to see if the UTR is fully contained within the coordinates of one exon
				singleexonUTR = False
				for exoncoord in exoncoords:
					exonstart, exonend = exoncoord[0], exoncoord[1]
					if exonstart <= UTRcoords[0] and exonend >= UTRcoords[1] and len(UTRcoords) > 0:
						singleexonUTR = True
						txUTRcoords[transcriptID] = [UTRcoords]
						if geneID not in gene2tx:
							gene2tx[geneID] = []
						gene2tx[geneID].append(transcriptID)
						tx2gene[transcriptID] = geneID

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

					txUTRcoords[transcriptID] = UTRexoncoords
					if geneID not in gene2tx:
						gene2tx[geneID] = []
					gene2tx[geneID].append(transcriptID)
					tx2gene[transcriptID] = geneID

		

	print 'Looked through {0} genes for {1} transcripts. Found 3\' UTRs for {2} of them.'.format(genecount, len(txs), len(txUTRcoords))


	#Output to gff
	with open(outputgff, 'w') as f:
		#txUTRcoords = {} #{ENSMUST : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
		#txchrmstrand = {} #{ENSMUST: [chrm, strand]}
		#genechrmstrand = {} #{ENSMUSG : [chrm , strand]}
		#gene2tx = {} # {ENSMUSG : ENSMUST}
		for gene in gene2tx:

			#Check to see if this gene has at least one transcript with a 3' UTR
			has3UTR = False
			for transcript in gene2tx[gene]:
				if len(txUTRcoords[transcript]) > 0:
					has3UTR = True

			if not has3UTR:
				continue

			genechrm = genechrmstrand[gene][0]
			genestrand = genechrmstrand[gene][1]
			genestart = geneboundaries[gene][0]
			genestop = geneboundaries[gene][1]
			if gene in e2sdict:
				geneshortname = e2sdict[gene]
			else:
				geneshortname = gene

			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(gene, geneshortname, gene)
			f.write(('\t').join([genechrm, 'gene', 'gene', str(genestart), str(genestop), '.', genestrand, '.', IDline]) + '\n')

			for transcript in gene2tx[gene]:
				#If it doesn't have a UTR, skip it
				if len(txUTRcoords[transcript]) == 0:
					continue
				exoncounter = 0
				txchrm = txchrmstrand[transcript][0]
				txstrand = txchrmstrand[transcript][1]
				UTRstart = txUTRcoords[transcript][0][0]
				UTRstop = txUTRcoords[transcript][-1][1]
				IDline = 'ID=UTR3:{0};Parent=gene:{1};gene_id={2}'.format(transcript, gene, gene)
				f.write(('\t').join([txchrm, 'UTR3', 'UTR3', str(UTRstart), str(UTRstop), '.', txstrand, '.', IDline]) + '\n')
				for UTRexon in txUTRcoords[transcript]:
					exoncounter +=1
					exonstart = UTRexon[0]
					exonstop = UTRexon[1]
					IDline = 'ID=exon:{0}.utr3exon{1};Parent=UTR3:{2}'.format(transcript, exoncounter, transcript)
					f.write(('\t').join([txchrm, 'exon', 'exon', str(exonstart), str(exonstop), '.', txstrand, '.', IDline]) + '\n')

		'''
		for transcript in txUTRcoords:
			#If there is no UTR here, skip it
			if len(txUTRcoords[transcript]) == 0:
				continue
			exoncounter = 0
			transcriptID = transcript.split('_')[0]
			chrm = transcript.split('_')[1]
			strand = transcript.split('_')[2]
			#transcriptID = ('_').join([transcript.split('_')[0], transcript.split('_')[1]])
			#chrm = transcript.split('_')[2]
			#strand = transcript.split('_')[3]
			geneID = tx2gene[transcriptID]
			if geneID in e2sdict:
				geneshortname = e2sdict[geneID]
			else:
				geneshortname = geneID
			genestart, geneend = geneboundaries[geneID][0], geneboundaries[geneID][1]
			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(geneID, geneshortname, geneID)
			f.write(('\t').join([chrm, 'gene', 'gene', str(genestart), str(geneend), '.', strand, '.', IDline]) + '\n')
			UTRstart = txUTRcoords[transcript][0][0]
			UTRstop = txUTRcoords[transcript][-1][1]
			IDline = 'ID=UTR3:{0};Parent=gene:{1};gene_id={2}'.format(transcriptID, geneID, geneID)
			f.write(('\t').join([chrm, 'UTR3', 'UTR3', str(UTRstart), str(UTRstop), '.', strand, '.', IDline]) + '\n')
			for UTRexon in txUTRcoords[transcript]:
				exoncounter +=1
				exonstart = UTRexon[0]
				exonstop = UTRexon[1]
				IDline = 'ID=exon:{0}.utr3exon{1};Parent=UTR3:{2}'.format(transcriptID, exoncounter, transcriptID)
				f.write(('\t').join([chrm, 'exon', 'exon', str(exonstart), str(exonstop), '.', strand, '.', IDline]) + '\n')
		'''

	return txUTRcoords, tx2gene, txchrmstrand


def get5UTRcoords(gff, ens2short, txs, outputgff):
	txUTRcoords = {} #{ENSMUST_chrm_strand : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
	txchrmstrand = {} #{ENSMUST: [chrm, strand]}
	genechrmstrand = {} #{ENSMUSG : [chrm , strand]}
	gene2tx = {} # {ENSMUSG : [ENSMUST]}
	tx2gene = {} # {ENSMUST : ENSMSUG}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
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

	a = []
	for tx in txs:
		a.append(tx.split('.')[0])
	txs = a

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1
		if genecount % 10000 == 0:
			print 'Gene {0}...'.format(genecount)
		geneID = str(gene.id).replace('gene:', '')
		geneboundaries[geneID] = [gene.start, gene.end]
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start'):
			#If this transcript has no coding exons, skip it:
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) == 0:
				continue

			transcriptID = str(transcript.id).replace('transcipt:', '').split('.')[0]
			if transcriptID in txs:

				if geneID not in genechrmstrand:
					genechrmstrand[geneID] = [gene.chrom, gene.strand]

				if transcriptID not in txchrmstrand:
					txchrmstrand[transcriptID] = [transcript.chrom, transcript.strand]


				tx2gene[transcriptID] = geneID
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
					#If the transcript starts right where the CDS starts, then there's no UTR
					if CDSstart == transcript.start:
						continue
					UTR5end = CDSstart - 1
					UTRcoords = [transcript.start, UTR5end]
				elif transcript.strand == '-':
					CDSstart = max(CDScoords, key = itemgetter(1))[1]
					#If the transcript starts right where the CDS starts, then there's no UTR
					if CDSstart == transcript.end:
						continue
					UTR5end = CDSstart + 1
					UTRcoords = [UTR5end, transcript.end]

				#Check to see if the UTR is fully contained within the coordinates of one exon
				singleexonUTR = False
				for exoncoord in exoncoords:
					exonstart, exonend = exoncoord[0], exoncoord[1]
					if exonstart <= UTRcoords[0] and exonend >= UTRcoords[1]:
						singleexonUTR = True
						txUTRcoords[transcriptID] = [UTRcoords]
						if geneID not in gene2tx:
							gene2tx[geneID] = []
						gene2tx[geneID].append(transcriptID)
						tx2gene[transcriptID] = geneID

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

					txUTRcoords[transcriptID] = UTRexoncoords
					if geneID not in gene2tx:
						gene2tx[geneID] = []
					gene2tx[geneID].append(transcriptID)
					tx2gene[transcriptID] = geneID

	print 'Looked through {0} genes for {1} transcripts. Found 5\' UTRs for {2} of them.'.format(genecount, len(txs), len(txUTRcoords))


	#Output to gff
	#txUTRcoords = {} #{ENSMUST : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
	#txchrmstrand = {} #{ENSMUST: [chrm, strand]}
	#genechrmstrand = {} #{ENSMUSG : [chrm , strand]}
	#gene2tx = {} # {ENSMUSG : ENSMUST}
	with open(outputgff, 'w') as f:
		for gene in gene2tx:

			#Check to see if this gene has at least one transcript with a 3' UTR
			has5UTR = False
			for transcript in gene2tx[gene]:
				if len(txUTRcoords[transcript]) > 0:
					has5UTR = True

			if not has5UTR:
				continue

			genechrm = genechrmstrand[gene][0]
			genestrand = genechrmstrand[gene][1]
			genestart = geneboundaries[gene][0]
			genestop = geneboundaries[gene][1]
			if gene in e2sdict:
				geneshortname = e2sdict[gene]
			else:
				geneshortname = gene

			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(gene, geneshortname, gene)
			f.write(('\t').join([genechrm, 'gene', 'gene', str(genestart), str(genestop), '.', genestrand, '.', IDline]) + '\n')

			for transcript in gene2tx[gene]:
				#If it doesn't have a UTR, skip it
				if len(txUTRcoords[transcript]) == 0:
					continue
				exoncounter = 0
				txchrm = txchrmstrand[transcript][0]
				txstrand = txchrmstrand[transcript][1]
				UTRstart = txUTRcoords[transcript][0][0]
				UTRstop = txUTRcoords[transcript][-1][1]
				IDline = 'ID=UTR5:{0};Parent=gene:{1};gene_id={2}'.format(transcript, gene, gene)
				f.write(('\t').join([txchrm, 'UTR5', 'UTR5', str(UTRstart), str(UTRstop), '.', txstrand, '.', IDline]) + '\n')
				for UTRexon in txUTRcoords[transcript]:
					exoncounter +=1
					exonstart = UTRexon[0]
					exonstop = UTRexon[1]
					IDline = 'ID=exon:{0}.utr5exon{1};Parent=UTR5:{2}'.format(transcript, exoncounter, transcript)
					f.write(('\t').join([txchrm, 'exon', 'exon', str(exonstart), str(exonstop), '.', txstrand, '.', IDline]) + '\n')


	'''
	with open(outputgff, 'w') as f:
		#txUTRcoords = {} #{ENSMUST_chrm_strand : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
		for transcript in txUTRcoords:
			exoncounter = 0
			#transcriptID = transcript.split('_')[0]
			#chrm = transcript.split('_')[1]
			#strand = transcript.split('_')[2]
			geneID = tx2gene[transcriptID]
			geneshortname = e2sdict[geneID]
			genestart, geneend = geneboundaries[geneID][0], geneboundaries[geneID][1]
			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(geneID, geneshortname, geneID)
			f.write(('\t').join([chrm, 'gene', 'gene', str(genestart), str(geneend), '.', strand, '.', IDline]) + '\n')
			UTRstart = txUTRcoords[transcript][0][0]
			UTRstop = txUTRcoords[transcript][-1][1]
			IDline = 'ID=UTR5:{0};Parent={1};gene_id={2}'.format(transcriptID, geneID, geneID)
			f.write(('\t').join([chrm, 'UTR5', 'UTR5', str(UTRstart), str(UTRstop), '.', strand, '.', IDline]) + '\n')
			for UTRexon in txUTRcoords[transcript]:
				exoncounter +=1
				exonstart = UTRexon[0]
				exonstop = UTRexon[1]
				IDline = 'ID=exon:{0}.utr5exon{1};Parent=UTR5:{2}'.format(transcriptID, exoncounter, transcriptID)
				f.write(('\t').join([chrm, 'exon', 'exon', str(exonstart), str(exonstop), '.', strand, '.', IDline]) + '\n')

	'''

	return txUTRcoords, tx2gene, txchrmstrand

def getSequences(regioncoords, genomefasta, ens2short, tx2gene, txchrmstrand, outfasta):
	#txUTRcoords = {} #{ENSMUST : [[utrexon1start, utrexon1stop], [utrexon2start, utrexon2stop]]}
	#txchrmstrand = {} #{ENSMUST: [chrm, strand]}
	#tx2gene = {} #{ENSMUST : ENSMUSG}
	#ens2short for mm10 is Ensembl_to_genename.txt
	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	seqs = {} #{txname : CDSseq}
	e2sdict = {} #{ENSGene : shortname}
	chrmswithoutseq = [] #Chromsome names that are in coords but that don't have a fasta entry in genomefasta

	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[2]
	infh.close()

	for tx in regioncoords:
		seq = ''
		txname = tx
		chrm = txchrmstrand[tx][0]
		strand = txchrmstrand[tx][1]
		#txname = ('_').join([tx.split('_')[0], tx.split('_')[1]])
		#chrm = tx.split('_')[2]
		#strand = tx.split('_')[3]

		#Is this chromosome in genomefasta?
		if chrm not in seq_dict: 
			if chrm not in chrmswithoutseq:
				print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
				chrmswithoutseq.append(chrm)
			continue

		for coords in regioncoords[tx]:
			start = coords[0]
			end = coords[1]
			if strand == '+':
				exonseq = seq_dict[chrm].seq[start-1:end].upper()
				seq += exonseq
			elif strand == '-':
				exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
				newseq = exonseq + seq
				seq = newseq

		genename = tx2gene[txname].split('.')[0]
		if genename in e2sdict:
			shortname = e2sdict[genename]
			seqs[txname + '_' + genename + '_' + shortname] = str(seq)
		else:
			shortname = genename
			seqs[txname + '_' + genename + '_' + shortname] = str(seq)

	with open(outfasta, 'w') as f:
		for seq in seqs:
			f.write('>' + seq + '\n' + str(seqs[seq]) + '\n')






if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--transcripts', type = str, help = 'List of ensembl transcript IDs.')
	parser.add_argument('--gff', type = str, help = 'Genome annotation containing transcript IDs.')
	parser.add_argument('--ens2short', type = str, help = 'File containing ensembl IDs to gene short name relations. Usually Ensembl_to_genename.txt')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--outputgff', type = str, help = 'Output file of gene regions in gff format.')
	parser.add_argument('--outputfasta', type = str, help = 'Output file of gene regions in fasta format.')
	parser.add_argument('--region', type = str, choices = ['UTR5', 'CDS', 'UTR3'])
	args = parser.parse_args()

	txs = []
	with open(args.transcripts, 'r') as f:
		for line in f:
			line = line.strip()
			txs.append(line)

	if args.region == 'CDS':
		coords, tx2gene, txchrmstrand = getCDScoords(args.gff, args.ens2short, txs, args.outputgff)
		getSequences(coords, args.genomefasta, args.ens2short, tx2gene, txchrmstrand, args.outputfasta)

	elif args.region == 'UTR3':
		coords, tx2gene, txchrmstrand = get3UTRcoords(args.gff, args.ens2short, txs, args.outputgff)
		getSequences(coords, args.genomefasta, args.ens2short, tx2gene, txchrmstrand, args.outputfasta)

	elif args.region == 'UTR5':
		coords, tx2gene, txchrmstrand = get5UTRcoords(args.gff, args.ens2short, txs, args.outputgff)
		getSequences(coords, args.genomefasta, args.ens2short, tx2gene, txchrmstrand, args.outputfasta)





