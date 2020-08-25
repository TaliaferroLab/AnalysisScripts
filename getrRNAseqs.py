import gffutils
from Bio import SeqIO
import gzip
import argparse
import os

def getrRNAcoords(gff):
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	txs = db.features_of_type('transcript')

	txcount = 0
	rRNAcount = 0
	rRNAcoords = {} #{txid : [[exon1start, exon1stop], [exon2start, exon2stop]]}
	for tx in txs:
		txcount +=1
		if txcount % 10000 == 0:
			print 'Transcript {0}...'.format(txcount)
		biotype = tx.attributes['gene_type'][0]
		if biotype == 'rRNA':
			rRNAcount +=1
			txid = str(tx.id) + '_' + tx.chrom + '_' + tx.strand
			rRNAcoords[txid] = []
			for exon in db.children(tx, featuretype = 'exon', level = 1):
				rRNAcoords[txid].append([exon.start, exon.end])

	print 'Looked through {0} transcripts and found {1} rRNA transcripts.'.format(txcount, rRNAcount)

	return rRNAcoords

def getrRNAseqs(rRNAcoords, genomefasta, outputfasta):
	outfh = open(outputfasta, 'w')
	outfh.close()

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	chrmswithoutseq = [] #Chromsome names that are in rRNAcoords but that don't have a fasta entry in genomefasta

	for tx in rRNAcoords:
		txid = tx.split('_')[0]
		chrm = tx.split('_')[1]
		strand = tx.split('_')[2]

		rRNAseq = ''
		#Is this chromosome in genomefasta?
		if chrm not in seq_dict: 
			if chrm not in chrmswithoutseq:
				print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
				chrmswithoutseq.append(chrm)
			continue

		if strand == '+':
			for exon in rRNAcoords[tx]:
				exonstart = exon[0]
				exonend = exon[1]
				exonseq = seq_dict[chrm].seq[exonstart-1:exonend].upper()
				rRNAseq += exonseq

		elif strand == '-':
			for exon in rRNAcoords[tx]:
				exonstart = exon[0]
				exonend = exon[1]
				exonseq = seq_dict[chrm].seq[exonstart-1:exonend].reverse_complement().upper()
				newseq = exonseq + rRNAseq
				rRNAseq = newseq

		with open(outputfasta, 'a') as f:
			f.write('>' + txid + '\n' + str(rRNAseq) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mm10gff', type = str, help = 'mm10 annotation in gff format.  Usually gencodecomprehensive.mm10.gff3.')
	parser.add_argument('--mm10genome', type = str, help = 'mm10 genome in fasta format.  Usually mm10.fasta.gz')
	parser.add_argument('--outputfasta', type = str, help = 'rRNA sequences in fasta format.')
	args = parser.parse_args()

	rRNAcoords = getrRNAcoords(args.mm10gff)
	getrRNAseqs(rRNAcoords, args.mm10genome, args.outputfasta)




