#Given a fasta file, get the length and GC content for each seq
from Bio import SeqIO
import argparse

def getGCcontent(seq):
	seqlen = len(seq)
	GCcount = seq.count('G') + seq.count('C')
	GCcontent = GCcount / float(seqlen)

	return GCcontent

def iteratefasta(fasta):
	d = {} #{seqname : [len, gccontent]}
	for record in SeqIO.parse(fasta, 'fasta'):
		ID = str(record.id).split('.')[0]
		seqlen = len(str(record.seq))
		gccontent = getGCcontent(str(record.seq))
		d[ID] = [seqlen, gccontent]

	return d

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file.')
	parser.add_argument('--region', type = str, help = 'Transcript region.')
	parser.add_argument('--output', type = str, help = 'Output file.')
	args = parser.parse_args()

	with open(args.output, 'w') as outfh:
		outfh.write(('\t').join(['Gene', 'region', 'length', 'gccontent']) + '\n')
		d = iteratefasta(args.fasta)
		for seq in d:
			seqlen = d[seq][0]
			gc = d[seq][1]
			outfh.write(('\t').join([seq.split('.')[0], args.region, str(seqlen), str(gc)]) + '\n')
