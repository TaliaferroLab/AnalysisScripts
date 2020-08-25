#For quantification with ribomap, we need to make a fasta of all transcripts that have a valid CDS.
#I'm going to define having a valid CDS has having start and stop codons. Some transcripts have
#features marked as "CDS", but when you stitch exons together, you find that they don't have the
#proper terminal codons.  

#In addition to a fasta file, we will also need a file telling where UTR/CDS boundaries are for each transcript.

import gffutils
import os
from Bio import SeqIO
import gzip
import sys
import argparse

#First, make all CDS sequences to figure out which ones are valid

def getCDScoords(gff):
	CDScoords = {} #{ENSTID_chrm_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	genecounter = 0
	txcounter = 0
	geneswithcdscounter = 0
	txwithcdscounter = 0
	for gene in genes:
		genecounter +=1
		if genecounter % 10000 == 0:
			print 'Gene {0}...'.format(genecounter)
		genehasCDS = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'start'):
			txcounter +=1
			if len(list(db.children(transcript, featuretype = 'CDS'))) == 0:
				continue
			genehasCDS = True
			txwithcdscounter +=1
			txid = str(transcript.id)
			chrm = str(transcript.chrom)
			strand = str(transcript.strand)
			ID = txid + '_' + chrm + '_' + strand
			CDScoords[ID] = []
			for cdsexon in db.children(transcript, featuretype = 'CDS', level = 1, order_by = 'start'):
				exoncoords = [cdsexon.start, cdsexon.end]
				CDScoords[ID].append(exoncoords)

		if genehasCDS:
			geneswithcdscounter +=1

	print 'Looked through {0} genes. {1} had coding exons.'.format(genecounter, geneswithcdscounter)
	print 'Looked through {1} transcripts. {1} had coding exons.'.format(txcounter, txwithcdscounter)

	return CDScoords


#Now figure out which ones have valid terminal codons
def getSequences(CDScoords, genomefasta):
	chrmswithoutseq = [] #Chromsome names that are in allCDScoords but that don't have a fasta entry in genomefasta
	validcdstxs = [] #list of tx names that have a valid cds
	txcounter = 0
	seqcounter = 0
	tooshort = 0
	wrongends = 0


	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	
	for tx in CDScoords:
		txcounter +=1
		CDSseq = ''
		txname = tx.split('_')[0]
		chrm = tx.split('_')[1]
		strand = tx.split('_')[2]

		#Is this chromosome in genomefasta?
		if chrm not in seq_dict:
			if chrm not in chrmswithoutseq:
				print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
				chrmswithoutseq.append(chrm)
			continue

		for coords in CDScoords[tx]:
			start = coords[0]
			end = coords[1]
			if strand == '+':
				exonseq = seq_dict[chrm].seq[start-1:end].upper()
				CDSseq += exonseq
			elif strand == '-':
				exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
				newseq = exonseq + CDSseq
				CDSseq = newseq

		if len(CDSseq) < 20:
			tooshort +=1
			continue

		if CDSseq[0:3] == 'ATG' and (CDSseq[-3:] == 'TAA' or CDSseq[-3:] == 'TGA' or CDSseq[-3:] == 'TAG'):
			seqcounter +=1
			validcdstxs.append(txname)
				
		else:
			wrongends +=1
			continue

	print 'Had CDS seqs for {0} transcripts. {1} were valid CDS sequences.  {2} were too short and {3} had incorrect terminal codons.'.format(txcounter, seqcounter, tooshort, wrongends)

	return validcdstxs

#Write the sequences of transcripts that have valid CDS sequences
def makefasta(gff, validcdstxs, genomefasta, outputfasta):
	txcoords = {} #{ENSTID_chrm_strand : [[exon1start, exon1end], [exon2start, exon2end]]}
	genecounter = 0
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecounter +=1
		if genecounter % 10000 == 0:
			print 'Gene {0}...'.format(genecounter)
		chrm = gene.chrom
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if str(transcript.id) in validcdstxs:
				txname = str(transcript.id).split('_')[0].replace('transcript:', '').split('.')[0]
				ID = txname + '_' + chrm + '_' + strand
				exoncoords = []
				for exon in db.children(transcript, featuretype = 'exon', level = 1, order_by = 'start'):
					coords = [exon.start, exon.end]
					exoncoords.append(coords)
				txcoords[ID] = exoncoords

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	chrmswithoutseq = [] #Chromsome names that are in allCDScoords but that don't have a fasta entry in genomefasta
	with open(outputfasta, 'w') as outfh:
		for tx in txcoords:
			txname = tx.split('_')[0]
			chrm = tx.split('_')[1]
			strand = tx.split('_')[2]
			seq = ''

			#Is this chromosome in genomefasta?
			if chrm not in seq_dict:
				if chrm not in chrmswithoutseq:
					print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
					chrmswithoutseq.append(chrm)
				continue

			for coords in txcoords[tx]:
				start = coords[0]
				end = coords[1]
				if strand == '+':
					exonseq = seq_dict[chrm].seq[start-1:end].upper()
					seq += exonseq
				elif strand == '-':
					exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
					newseq = exonseq + seq
					seq = newseq

			outfh.write('>' + txname + '\n' + str(seq) + '\n')

#Make a gff containing only genes with at least one valid transcript
#and all valid transcripts

def makegff(gff, validcdstxs, outputgff):
	print 'Making gff...'
	outfh = open(outputgff, 'w')
	outfh.close()

	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genehasvalidCDS = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if str(transcript.id) in validcdstxs:
				genehasvalidCDS = True
				with open(outputgff, 'a') as outfh:
					print >> outfh, gene
				break
		
		if not genehasvalidCDS:
			continue
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if str(transcript.id) in validcdstxs:
				with open(outputgff, 'a') as outfh:
					print >> outfh, transcript

					for feature in db.children(transcript):
						print >> outfh, feature


#Get the fasta file CDS/UTR boundaries
def getCDSboundaries(gff, validcdstxs, txfasta, outfile):
	txfastadict = {}
	for record in SeqIO.parse(txfasta, 'fasta'):
		seq = str(record.seq)
		txfastadict[str(record.id)] = seq

	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')
	for gene in genes:
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if str(transcript.id) in validcdstxs:
				exoncoords = []
				cdscoords = []
				if transcript.strand == '+':
					for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
						exonrange = range(exon.start, exon.end + 1)
						exoncoords += exonrange
					for cds in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
						cdsrange = range(cds.start, cds.end + 1)
						cdscoords += cdsrange
				elif transcript.strand == '-':
					for exon in db.children(transcript, featuretype = 'exon', order_by = 'start', reverse = True):
						exonrange = range(exon.start, exon.end + 1)
						exonrange.reverse()
						exoncoords += exonrange
					for cds in db.children(transcript, featuretype = 'CDS', order_by = 'start', reverse = True):
						cdsrange = range(cds.start, cds.end + 1)
						cdsrange.reverse()
						cdscoords += cdsrange

				cdsstart = cdscoords[0]
				cdsend = cdscoords[-1]
				#First nt of start codon
				startindex = exoncoords.index(cdsstart)
				#One after last nt of stop codon
				stopindex = exoncoords.index(cdsend) + 1
				txshortname = str(transcript.id).split('.')[0]
				startseq = txfastadict[txshortname][startindex : startindex + 3]
				stopseq = txfastadict[txshortname][stopindex- 3 : stopindex]
				correctindex = False
				if startseq == 'ATG' and (stopseq == 'TAA' or stopseq == 'TGA' or stopseq == 'TAG'):
					correctindex = True
					#print transcript.id, transcript.strand, cdsstart, cdsend, startindex, stopindex, startseq, stopseq, 'RIGHT'

					#This is the CDS range file
					with open(outfile, 'a') as outfh:
						outfh.write(txshortname + '\t' + str(startindex) + '\t' + str(stopindex) + '\n')
				if not correctindex:
					print transcript.id, transcript.strand, cdsstart, cdsend, startindex, stopindex, startseq, stopseq, 'WRONG'





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF genome annotation file. Usually gencodecomprehensive.mm10.gff3.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in gzipped fasta format.')
	parser.add_argument('--txfasta', type = str, help = 'Output fasta file of valid CDS-contatining transcripts.')
	parser.add_argument('--txgff', type = str, help = 'Output gff file of valid CDS-contatining transcripts.')
	parser.add_argument('--CDSrangefile', type = str, help = 'Output CDS range file.')
	args = parser.parse_args()

	CDScoords = getCDScoords(args.gff)
	validcdstxs = getSequences(CDScoords, args.genomefasta)
	makefasta(args.gff, validcdstxs, args.genomefasta, args.txfasta)
	makegff(args.gff, validcdstxs, args.txgff)
	getCDSboundaries(args.gff, validcdstxs, args.txfasta, args.CDSrangefile)


