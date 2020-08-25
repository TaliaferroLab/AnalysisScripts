#Given a list of genes and a genome annotation, go through and figure out which ones have 0 coding transcripts.
#Then, for those genes, get the longest annotated transcript.
#Report this transcript in both gff and fasta form.


import gffutils
from Bio import SeqIO
import operator
import sys
import os
import gzip
import argparse

def seeifcoding(gene, db):
	#See if this gene has any coding transcripts
	coding = False
	for transcript in db.children(gene, featuretype = 'transcript', level = 1):
		if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
			coding = True

	return coding

def getlongesttranscript(gene, db):
	txlengths = {} #{txid : length of exons}
	for transcript in db.children(gene, featuretype = 'transcript', level = 1):
		transcriptid = str(transcript.id)
		tlength = 0
		for exon in db.children(transcript, featuretype = 'exon', level = 1):
			elength = exon.end - exon.start
			tlength += elength

		txlengths[transcriptid] = tlength
		longesttx = max(txlengths.iteritems(), key = operator.itemgetter(1))[0]
	return longesttx

def getnoncodingtxs(gff, genesofint):
	#Given a list of genes, figure out which ones have 0 coding transcripts.
	#Then go through and get the longest transcript for that gene.
	geneandtx = {} #{ENSMUSG : ENSMUST}
	geneboundaries = {} #{ENSMUSG : [genestart, genestop]}
	genecount = 0
	noncodinggenecount = 0

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
		geneid = str(gene.id).replace('gene:', '').split('.')[0]
		if geneid in genesofint:
			genecount +=1
			iscoding = seeifcoding(gene, db)
			if not iscoding:
				noncodinggenecount +=1
				longesttx = getlongesttranscript(gene, db)
				geneandtx[geneid] = longesttx

	print 'Looked for {0} genes in gff and found {1} of them. Of these {2} had no coding transcripts.'.format(len(genesofint), genecount, noncodinggenecount)

	return geneandtx

def getcoords(geneandtx, gff, outputgff):
	geneboundaries = {} #{ENSMUSG : [chrm, genestart, genestop, strand]}
	txchrmstrand = {} #{txid : [chrm, txstart, txstop, strand]}
	txcoords = {} #{ENSMUSG : {ENSMUST : [[exon1start, exon1stop], [exon2start, exon2stop]]}}
	#Given a list of transcript, return exonic coords
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
		geneid = str(gene.id).replace('gene:', '').split('.')[0]
		if geneid in geneandtx:
			geneboundaries[geneid] = [gene.chrom, gene.start, gene.end, gene.strand]
			longesttx = geneandtx[geneid]
			txcoords[geneid] = {}
			for transcript in db.children(gene, featuretype = 'transcript', level = 1):
				if str(transcript.id) == longesttx:
					txchrmstrand[str(transcript.id)] = [transcript.chrom, transcript.start, transcript.end, transcript.strand]
					exoncoords = []
					for exon in db.children(transcript, featuretype = 'exon', level = 1, order_by = 'start'):
						exoncoords.append([exon.start, exon.end])
					txcoords[geneid][str(transcript.id)] = exoncoords

	with open(outputgff, 'w') as f:
		for gene in txcoords:
			genechrm = geneboundaries[gene][0]
			genestart = geneboundaries[gene][1]
			geneend = geneboundaries[gene][2]
			genestrand = geneboundaries[gene][3]
			IDline = 'ID=gene:{0};Name={1};gene_id={2}'.format(gene, gene, gene)
			f.write(('\t').join([genechrm, 'noncoding', 'gene', str(genestart), str(geneend), '.', genestrand, '.', IDline]) + '\n')
			tx = txcoords[gene].keys()[0]
			txchrm = txchrmstrand[tx][0]
			txstart = txchrmstrand[tx][1]
			txend = txchrmstrand[tx][2]
			txstrand = txchrmstrand[tx][3]
			IDline = 'ID=transcript:{0};Parent=gene:{1};gene_id={2}'.format(tx, gene, gene)
			f.write(('\t').join([txchrm, 'noncoding', 'transcript', str(txstart), str(txend), '.', txstrand, '.', IDline]) + '\n')
			exoncounter = 0
			for exon in txcoords[gene][tx]:
				exoncounter +=1
				exonstart = exon[0]
				exonend = exon[1]
				IDline = 'ID=exon:{0}.exon{1};Parent=transcript:{2};gene_id={3}'.format(tx, exoncounter, tx, gene)
				f.write(('\t').join([txchrm, 'noncoding', 'exon', str(exonstart), str(exonend), '.', txstrand, '.', IDline]) + '\n')

	return txcoords, txchrmstrand

def getSeqs(txcoords, txchrmstrand, genomefasta, outputfasta):
	#txcoords = {} #{ENSMUSG : {ENSMUST : [[exon1start, exon1stop], [exon2start, exon2stop]]}}
	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	seqs = {} #{genename : {txname : cDNAseq}}

	for gene in txcoords:
		seqs[gene] = {}
		for tx in txcoords[gene]:
			seq = ''
			txname = txcoords[gene].keys()[0]
			chrm = txchrmstrand[txname][0]
			strand = txchrmstrand[txname][3]
			for exon in txcoords[gene][tx]:
				start = exon[0]
				end = exon[1]
				if strand == '+':
					exonseq = seq_dict[chrm].seq[start-1:end].upper()
					seq += exonseq
				elif strand == '-':
					exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
					seq = exonseq + seq

			seqs[gene][txname] = seq

	with open(outputfasta, 'w') as f:
		for gene in seqs:
			for tx in seqs[gene]:
				seq = seqs[gene][tx]
				f.write('>' + gene + '_' + tx + '\n' + str(seq) + '\n')



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--genes', type = str, help = 'List of ENSEMBL gene IDs to look through.')
	parser.add_argument('--gff', type = str, help = 'Genome annotation.')
	parser.add_argument('--outputgff', type = str, help = 'Output gff of the longest transcript of any noncoding gene in the list.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in gzipped fasta format.')
	parser.add_argument('--outputfasta', type = str, help = 'Output fasta file of longest transcript of any noncoding gene in the list.')
	args = parser.parse_args()


	genesofint = []
	with open(args.genes, 'r') as f:
		for line in f:
			line = line.strip()
			genesofint.append(line)

	geneandtx = getnoncodingtxs(args.gff, genesofint)
	txcoords, txchrmstrand = getcoords(geneandtx, args.gff, args.outputgff)
	getSeqs(txcoords, txchrmstrand, args.genomefasta, args.outputfasta)



