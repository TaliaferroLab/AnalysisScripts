import gffutils
import os
import sys
import argparse
from Bio import SeqIO
import gzip

#python3

def makeregions(gff, genomefasta, outfile):
	upstreamantiseqs = {} #{genename : seq}
	downstreampolyAseqs = {} #{genename : seq}

	print('Indexing genome sequence...')
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta, 'rt'), 'fasta'))
	print('Done indexing!')

	#Make gff databases
	print('Indexing gff...')
	gff_fn = os.path.abspath(gff)
	db_fn = gff_fn + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)
	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')

	genecounter = 0
	with open(outfile, 'w') as outfh:
		for gene in genes:
			if 'protein_coding' not in gene.attributes['gene_type']:
				continue
			genecounter +=1
			if genecounter % 1000 == 0:
				print('Gene {0}...'.format(genecounter))

			outfh.write(str(gene) + '\n')

			if gene.strand == '+':
				#upstream antisense
				start = gene.start - 200
				end = gene.start - 1
				outfh.write(('\t').join([gene.chrom, 'gencode28', 'upstreamantisense', str(start), str(end), '.', '-', '.', 'ID={0}.upstreamantisense;gene_id={0};Parent={0}'.format(gene.id)]) + '\n')
				upstreamantiseq = seq_dict[gene.chrom].seq[start : end].reverse_complement().upper()

				#downstream of polyA
				start = gene.end + 100
				end = gene.end + 5000
				outfh.write(('\t').join([gene.chrom, 'gencode28', 'downstreampolyA', str(start), str(end), '.', '+', '.', 'ID={0}.downstreampolyA;gene_id={0};Parent={0}'.format(gene.id)]) + '\n')
				downstreampolyAseq = seq_dict[gene.chrom].seq[start : end].upper()

			elif gene.strand == '-':
				#upstream antisense
				start = gene.end + 1
				end = gene.end + 200
				outfh.write(('\t').join([gene.chrom, 'gencode28', 'upstreamantisense', str(start), str(end), '.', '+', '.', 'ID={0}.upstreamantisense;gene_id={0};Parent={0}'.format(gene.id)]) + '\n')
				upstreamantiseq = seq_dict[gene.chrom].seq[start : end].upper()

				#downstream of polyA
				start = gene.start - 5000
				end = gene.start - 100
				outfh.write(('\t').join([gene.chrom, 'gencode28', 'downstreampolyA', str(start), str(end), '.', '-', '.', 'ID={0}.downstreampolyA;gene_id={0};Parent={0}'.format(gene.id)]) + '\n')
				downstreampolyAseq = seq_dict[gene.chrom].seq[start : end].reverse_complement().upper()

			upstreamantiseqs[str(gene.id).split('.')[0] + '_upstreamanti'] = upstreamantiseq
			downstreampolyAseqs[str(gene.id).split('.')[0] + '_downstreampolyA'] = downstreampolyAseq

	#Write seqs
	with open('upstreamantiseqs.fa', 'w') as outfh:
		for gene in upstreamantiseqs:
			outfh.write('>' + gene + '\n' + str(upstreamantiseqs[gene]) + '\n')

	with open('downstreampolyAseqs.fa', 'w') as outfh:
		for gene in downstreampolyAseqs:
			outfh.write('>' + gene + '\n' + str(downstreampolyAseqs[gene]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF annotation to use.')
	parser.add_argument('--outfile', type = str, help = 'GFF output file.')
	parser.add_argument('--genomefasta', type = str, help = 'Fasta of genome sequence.')
	args = parser.parse_args()

	makeregions(args.gff, args.genomefasta, args.outfile)


