#python3

import os
from Bio import SeqIO
import sys
import pandas as pd
from collections import defaultdict

def getoligos(oligoresults):
	#Take a results file produced by Fractionation.Rmd (CADoligoresults.txt, for example)
	#and extract oligo IDs based on different properties (neurite loc, soma loc, etc.)
	df = pd.read_csv(oligoresults, sep = '\t', header = 0, index_col = 'oligo')
	#Filter for sig neurite loc oligos
	neuriteloc = df.query('neuritelog2FC > 0 & padj < 0.1').index.tolist()
	#All other oligos that are not neurite loc
	nonneuriteloc = df.query('neuritelog2FC < 0 | padj >= 0.1').index.tolist()
	#Filter for sig neurite soma oligos
	somaloc = df.query('neuritelog2FC < 0 & padj < 0.1').index.tolist()
	#All other oligos that are not soma loc
	nonsomaloc = df.query('neuritelog2FC > 0 | padj >= 0.1').index.tolist()

	return neuriteloc, nonneuriteloc, somaloc, nonsomaloc


def getGCcontents(fasta):
	#Get oligo sequences
	with open(fasta, 'r') as infh, open('GCcontents.txt', 'w') as outfh:
		outfh.write('oligo' + '\t' + 'GCcontent' + '\n')
		for record in SeqIO.parse(infh, 'fasta'):
			oligoid = record.id
			#Remove first 20 and last 20 nt (PCR handles)
			seq = record.seq[20:-20]
			GCcontent = (seq.count('C') + seq.count('G')) / len(seq)
			outfh.write(oligoid + '\t' + str(GCcontent) + '\n')

def getnuccontents(fasta):
	with open(fasta, 'r') as infh, open('Nuccontents.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'Acontent', 'Gcontent', 'Ccontent', 'Ucontent']) + '\n')
		for record in SeqIO.parse(infh, 'fasta'):
			oligoid = record.id
			#Remove first 20 and last 20 nt (PCR handles)
			seq = record.seq[20:-20]
			Acontent = seq.count('A') / float(len(seq))
			Gcontent = seq.count('G') / float(len(seq))
			Ccontent = seq.count('C') / float(len(seq))
			Ucontent = seq.count('T') / float(len(seq))
			outfh.write(('\t').join([oligoid, str(Acontent), str(Gcontent), str(Ccontent), str(Ucontent)]) + '\n')

#GC contents of the longest 5' UTR, CDS, and 3' UTR of all genes
#Longest fastas were produced by getLongestGeneRegions.py using GenocdeM17 comprehensive gff
def getgeneGCcontents(longest5utrfasta, longestcdsfasta, longest3utrfasta):
	#Get genes that are present in all three fastas
	utr5seqs = {} #all genes in these files, {gene : utr5seq}
	cdsseqs = {}
	utr3seqs = {}

	genesinallfiles = []
	
	with open(longest5utrfasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			genename = record.id.split('.')[0]
			utr5seqs[genename] = str(record.seq)
	with open(longestcdsfasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			genename = record.id.split('.')[0]
			cdsseqs[genename] = str(record.seq)
	with open(longest3utrfasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			genename = record.id.split('.')[0]
			utr3seqs[genename] = str(record.seq)

	utr5genes = set(utr5seqs.keys())
	cdsgenes = set(cdsseqs.keys())
	utr3genes = set(utr3seqs.keys())

	for gene in utr5genes:
		if gene in cdsgenes and gene in utr3genes:
			genesinallfiles.append(gene)

	with open('geneGCcontents.txt', 'w') as outfh:
		outfh.write(('\t').join(['ensembl_gene_id', 'utr5gc', 'cdsgc', 'utr3gc']) + '\n')
		for gene in genesinallfiles:
			utr5seq = utr5seqs[gene]
			cdsseq = cdsseqs[gene]
			utr3seq = utr3seqs[gene]

			utr5gc = (utr5seq.count('G') + utr5seq.count('C')) / float(len(utr5seq))
			cdsgc = (cdsseq.count('G') + cdsseq.count('C')) / float(len(cdsseq))
			utr3gc = (utr3seq.count('G') + utr3seq.count('C')) / float(len(utr3seq))

			outfh.write(('\t').join([gene, str(utr5gc), str(cdsgc), str(utr3gc)]) + '\n')




#neuriteloc, nonneuriteloc, somaloc, nonsomaloc = getoligos(sys.argv[1])
#getGCcontents(sys.argv[2])

#getnuccontents(sys.argv[1])

getgeneGCcontents(sys.argv[1], sys.argv[2], sys.argv[3])