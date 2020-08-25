import argparse
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
from numpy import mean

def nameconversion(ens_to_short):
	#Given a file that has ENSG IDs and their corresponding short names, make a dictionary.
	ens2short = {} # {ENSMUSG000000 : gene_short_name}
	infh = open(ens_to_short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			ens2short[line[0]] = line[1]
	infh.close()
	return ens2short


def getlongestUTRs(UTRgff, ens2short):
	#Given a gff of UTRs (usualy mm9_ensGene.3putrs.gff), get the longest UTR for that gene.
	UTRs = {} # {ENSGENEID : [chrm, start, stop, strand]}
	infh = open(UTRgff, 'r')
	for line in infh:
		line = line.strip().split('\t')
		chrm = line[0]
		start = int(line[3])
		stop = int(line[4])
		strand = line[6]
		gene = line[8].split(';')[-1]
		if gene in ens2short and 'random' not in chrm:
			#Don't deal with any gene that doesn't have a short name or is on chr_random
			gene_short_name = ens2short[gene]
		else:
			continue
		length = stop - start
		if gene_short_name not in UTRs:
			UTRs[gene_short_name] = [chrm, start, stop, strand]
		elif gene_short_name in UTRs:
			currentlength = UTRs[gene_short_name][2] - UTRs[gene_short_name][1]
			if length > currentlength:
				UTRs[gene_short_name] = [chrm, start, stop, strand]

	infh.close()
	print 'Have UTRs for {0} genes.'.format(len(UTRs))

	return UTRs

def getsequences(UTRs, genomefasta, genes):
	genesofinterest = []
	GCs = []
	lengths = []
	infh = open(genes, 'r')
	for line in infh:
		line = line.strip()
		genesofinterest.append(line)
	infh.close()

	seqs = {} # {genename : UTR_sequence}
	sys.stderr.write('Indexing genome sequence...\n')
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	sys.stderr.write('{0} chromosomes indexed.\n'.format(len(seq_dict)))
	for UTR in UTRs:
		chrm = UTRs[UTR][0]
		start = UTRs[UTR][1]
		stop = UTRs[UTR][2]
		strand = UTRs[UTR][3]
		if strand == '+':
			UTRseq = seq_dict[chrm].seq[start - 1 : stop].upper()
		elif strand == '-':
			UTRseq = seq_dict[chrm].seq[start - 1 : stop].upper().reverse_complement()

		if UTR in genesofinterest:
			seqs[UTR] = str(UTRseq)
			GCs.append(GC(str(UTRseq)))
			lengths.append(len(str(UTRseq)))


	print 'Started with {0} genes. Found UTR sequences for {1} of them. Their average GC content is {2}%.'.format(len(genesofinterest), len(seqs), mean(GCs))
	print 'Their average length is {0}.'.format(mean(lengths))

	outfh.close()
	return seqs

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--ens2short', type = str, help = 'File of tab delimited ENSGENEIDs and gene short names.')
	parser.add_argument('--genes', type = str, help = 'List of genes for which you want the 3\' UTRs.')
	parser.add_argument('--UTRgff', type = str, help = '3\'UTR coordinates in gff format. mm9_ensGene.3putrs.gff, for example.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--output', type = str, help = 'Output file.')
	args = parser.parse_args()

	ens2short = nameconversion(args.ens2short)
	UTRs = getlongestUTRs(args.UTRgff, ens2short)
	seqs = getsequences(UTRs, args.genomefasta, args.genes)
	outfh = open(args.output, 'w')
	for UTR in seqs:
		outfh.write('>' + UTR + '\n')
		outfh.write(seqs[UTR] + '\n')
	outfh.close()



