from Bio import SeqIO
import argparse
import numpy as np
from scipy import stats

def getgenes(genelist):
	#Read in a set of genes
	genes = []
	with open(genelist, 'r') as infh:
		for line in infh:
			line = line.strip()
			genes.append(line)

	return genes

def getgcmeta(seq, numberofbins):
	#Given a sequence and the number of bins to divide it into for the metagene,
	#get the gc content of each bin

	seq = str(seq).upper()
	seqlen = len(seq)
	binwidth = seqlen / float(numberofbins)
	gcdict = {}

	for i in range(numberofbins):
		binnumber = i + 1
		binstart = int(round((i / float(numberofbins)) * seqlen))
		binend = int(round(((i + 1) / float(numberofbins)) * seqlen))
		binseq = seq[binstart : binend]
		gccount = binseq.count('G') + binseq.count('C')
		gccontent = gccount / float(len(binseq))
		gcdict[binnumber] = gccontent

	return gcdict

def iteratefasta(fasta, genes, numberofbins):
	#Go through a fasta of all seqs.  If genename is one we are interested in, get gc bins.
	foundgenes = 0

	allseqs_gcdict = {}
	for i in range(numberofbins):
		allseqs_gcdict[i + 1] = []

	for record in SeqIO.parse(fasta, 'fasta'):
		genename = record.id.split('_')[0]
		if genename in genes:
			foundgenes +=1
			seq = str(record.seq)
			if len(seq) < numberofbins:
				continue
			gcdict = getgcmeta(seq, numberofbins)
			for binnumber in gcdict:
				allseqs_gcdict[binnumber].append(gcdict[binnumber])

	print 'Found sequences for {0} of {1} genes.'.format(foundgenes, len(genes))

	#Take mean gc content in each bin across all genes.
	meangcdict = {}
	for binnumber in allseqs_gcdict:
		meangc = np.mean(allseqs_gcdict[binnumber])
		sdgc = np.std(allseqs_gcdict[binnumber])
		segc = stats.sem(allseqs_gcdict[binnumber])
		meangcdict[binnumber] = [meangc, sdgc, segc]

	return meangcdict

def writeoutput(meangcdict, region, seqclass, outputfile):
	#Region is UTR5, CDS, or UTR3
	#seqclass is FBSminusrodown, FBSplusrounchanged, etc.
	with open(outputfile, 'a') as outfh:
		for binnumber in meangcdict:
			meangc = meangcdict[binnumber][0]
			sdgc = meangcdict[binnumber][1]
			segc = meangcdict[binnumber][2]
			outfh.write(('\t').join([str(binnumber), str(meangc), str(sdgc), str(segc), region, seqclass]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--genes', type = str, help = 'List of genes to consider.')
	parser.add_argument('--fasta', type = str, help = 'List of all sequences to consider, including those in genes of interest (and others).')
	parser.add_argument('--numberofbins', type = int, help = 'Number of bins to divide metagene into.')
	parser.add_argument('--region', type = str, help = 'What transcript region are the seqs in fasta from? UTR5, CDS, UTR3')
	parser.add_argument('--seqclass', type = str, help = 'What class of genes are these?  Affected, unaffected, etc.')
	parser.add_argument('--outputfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	genes = getgenes(args.genes)
	meangcdict = iteratefasta(args.fasta, genes, args.numberofbins)
	writeoutput(meangcdict, args.region, args.seqclass, args.outputfile)


