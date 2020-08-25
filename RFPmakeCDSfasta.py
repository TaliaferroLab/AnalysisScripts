#Given a genome annotation, make a fasta of only CDS sequences, excluding the first 50 nt.

import gffutils
import os
from Bio import SeqIO
import gzip
import argparse

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

def getSequences(CDScoords, genomefasta, outputfasta):
	chrmswithoutseq = [] #Chromsome names that are in allCDScoords but that don't have a fasta entry in genomefasta
	txcounter = 0
	seqcounter = 0
	tooshort = 0
	wrongends = 0


	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	
	with open(outputfasta, 'w') as outfh:
		for tx in CDScoords:
			txcounter +=1
			CDSseq = ''
			txname = tx.split('_')[0].replace('transcript:', '').split('.')[0]
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

			if len(CDSseq) < 200:
				tooshort +=1
				continue

			if CDSseq[0:3] == 'ATG' and (CDSseq[-3:] == 'TAA' or CDSseq[-3:] == 'TGA' or CDSseq[-3:] == 'TAG'):
				seqcounter +=1
				#Take off first 50 nt
				outfh.write('>' + txname + '\n' + str(CDSseq[50:]) + '\n')
				
			else:
				wrongends +=1
				continue

	print 'Had CDS seqs for {0} transcripts. Wrote seqs for {1} of them.  {2} were too short and {3} had incorrect terminal codons.'.format(txcounter, seqcounter, tooshort, wrongends)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Gff annotation of genome.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--outputfasta', type = str, help = 'Output fasta.')
	args = parser.parse_args()

	CDScoords = getCDScoords(args.gff)
	getSequences(CDScoords, args.genomefasta, args.outputfasta)





