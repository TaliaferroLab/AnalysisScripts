import os
import sys
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import numpy as np
import argparse

#python3

def getgenesswith2seqs(fasta):
	#Given a fasta of unique UTR seqs, get IDs of genes that have 2 and only 2 unique utr seqs
	utrcounts = defaultdict(int)
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			geneid = record.id.split('_')[0]
			utrcounts[geneid] +=1

	geneswith2seqs = []
	for gene in utrcounts:
		if utrcounts[gene] == 2:
			geneswith2seqs.append(gene)

	print('Looked through UTRs of {0} genes. {1} had only 2 unique UTR seqs.'.format(len(utrcounts), len(geneswith2seqs)))

	return geneswith2seqs

def getAandButrseqs(fasta, geneswith2seqs):
	#Take the genes that only have 2 unique UTR seqs, and get the "A" (proximal) and "B" (distal) sequences
	#This fasta is again the same fasta of all the unique utr seqs
	Aseqs = {} #{gene : AUTRseq}
	Bseqs = {} #{gene : BUTRseq}
	geneswith2seqs = set(geneswith2seqs) #faster membership tests
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			geneid = record.id.split('_')[0]
			if geneid in geneswith2seqs:
				utrseqnumber = record.id.split('_')[1]
				if utrseqnumber == 'uniqueUTR0':
					#Aseq
					Aseqs[geneid] = str(record.seq)
				elif utrseqnumber == 'uniqueUTR1':
					#Bseq
					Bseqs[geneid] = str(record.seq)

	return Aseqs, Bseqs

def classifygenes(LABRAToutput):
	#Given a table of LABRAT psis and pvalues, classify genes as either having a sig +deltapsi, sig -deltapsi, or nonsigpsi
	df = pd.read_csv(LABRAToutput, sep = '\t', header = 0)
	posdpsigenes = []
	negdpsigenes = []
	nonsiggenes = []
	genetypes = {} #{gene : genetype (ALE, TUTR, mixed)}
	for row in df.itertuples():
		gene = row.Gene
		genetype = row.genetype
		genetypes[gene] = genetype
		'''
		N2A2dpsi = row.N2ANeurite2 - row.N2ASoma2
		N2A3dpsi = row.N2ANeurite3 - row.N2ASoma3
		CAD1dpsi = row.CADNeurite1 - row.CADSoma1
		CAD2dpsi = row.CADNeurite2 - row.CADSoma2
		CAD3dpsi = row.CADNeurite3 - row.CADSoma3
		dpsi = np.mean([N2A2dpsi, N2A3dpsi, CAD1dpsi, CAD2dpsi, CAD3dpsi])
		'''
		somapsi = np.mean([row.GFPctrlSomaA, row.GFPctrlSomaB, row.GFPctrlSomaC, row.GFPhypoxiaSomaA, row.GFPhypoxiaSomaB, row.GFPhypoxiaSomaC])
		neuritepsi = np.mean([row.GFPctrlNeuriteA, row.GFPctrlNeuriteB, row.GFPctrlNeuriteC, row.GFPhypoxiaNeuriteA, row.GFPhypoxiaNeuriteB, row.GFPhypoxiaNeuriteC])
		dpsi = neuritepsi - somapsi
		qval = row.qval
		if qval < 0.2 and dpsi > 0:
			posdpsigenes.append(gene)
		elif qval < 0.2 and dpsi < 0:
			negdpsigenes.append(gene)
		elif qval >= 0.2:
			nonsiggenes.append(gene)

	return posdpsigenes, negdpsigenes, nonsiggenes, genetypes


def getGCs(Aseqs, Bseqs, posdpsigenes, negdpsigenes, nonsiggenes, genetypes):
		with open('dpsiGCs.MNGFP.txt', 'w') as outfh:
			outfh.write(('\t').join(['gene', 'category', 'AseqGC', 'BseqGC', 'genetype']) + '\n')
			for gene in Aseqs:
				Aseq = Aseqs[gene]
				Bseq = Bseqs[gene]
				Agc = (Aseq.count('G') + Aseq.count('C')) / len(Aseq)
				Bgc = (Bseq.count('G') + Bseq.count('C')) / len(Bseq)
				if gene in posdpsigenes:
					category = 'posdpsi'
					genetype = genetypes[gene]
				elif gene in negdpsigenes:
					category = 'negdpsi'
					genetype = genetypes[gene]
				elif gene in nonsiggenes:
					category = 'nonsigdpsi'
					genetype = genetypes[gene]
				else: #it's not in the dPSI table
					continue
				outfh.write(('\t').join([gene, category, str(Agc), str(Bgc), genetype]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta of unique UTR sequences.  Produced by UTRuniqueseqs_PF.py')
	parser.add_argument('--dPSItable', type = str, help = 'Table of LABRAT dPSI values. This is looking for ProteinCodingTxOnly/LABRATpsis.vanilla.3end.txt.pval')
	args = parser.parse_args()

	geneswith2seqs = getgenesswith2seqs(args.fasta)
	Aseqs, Bseqs = getAandButrseqs(args.fasta, geneswith2seqs)
	posdpsigenes, negdpsigenes, nonsiggenes, genetypes = classifygenes(args.dPSItable)
	getGCs(Aseqs, Bseqs, posdpsigenes, negdpsigenes, nonsiggenes, genetypes)
