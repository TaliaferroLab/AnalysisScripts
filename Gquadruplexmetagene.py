import re
from Bio import SeqIO
import os
import argparse

def getmatchpositions(regexp, seq):
	matches = [(m.start(0)) for m in re.finditer(regexp, seq)]
	return matches

def regexpmeta(fasta, bins):
	motif = r'(?=([AU]GG(.{0,7})[AU]GG(.{0,7})[AU]GG(.{0,7})[AU]GG))'
	bins = float(bins)
	motiffrac = {} #{binnumber : [number of times you find motif in this bin, number of times you looked]}
	for i in range(int(bins)):
		motiffrac[i + 1] = [0, 0]

	for record in SeqIO.parse(fasta, 'fasta'):
		seq = str(record.seq.transcribe())
		for i in range(len(seq)):
			pos = i + 1
			posbin = int(round((pos / float(len(seq))) * bins))
			if posbin == 0:
				posbin = 1
			motiffrac[posbin][1] +=1
		matches = getmatchpositions(motif, seq)
		for match in matches:
			pos = match + 1
			posbin = int(round((pos / float(len(seq))) * bins))
			if posbin == 0:
				posbin = 1
			motiffrac[posbin][0] +=1

	return motiffrac

def writebins(motiffrac, region, seqclass, outfile):
	if not os.path.isfile(outfile):
		with open(outfile, 'w') as f:
			f.write(('\t').join(['Bin', 'frac', 'region', 'class']) + '\n')
	with open(outfile, 'a') as f:
		for posbin in motiffrac:
			frac = motiffrac[posbin][0] / float(motiffrac[posbin][1])
			f.write(('\t').join([str(posbin), str(frac), region, seqclass]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to look through.')
	parser.add_argument('--bins', type = str, help = 'Number of bins for metagene.')
	parser.add_argument('--region', type = str, help = 'Transcript region of fasta.')
	parser.add_argument('--seqclass', type = str, help = 'Delta LR class of fasta.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	motiffrac = regexpmeta(args.fasta, args.bins)
	writebins(motiffrac, args.region, args.seqclass, args.outfile)

