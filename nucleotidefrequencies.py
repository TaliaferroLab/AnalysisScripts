#Usage: python nucleotidefrequencies.py <fasta file> <output file>
#Output is tab delimited frequencies of A, G, C, U

from Bio import SeqIO
import sys

def getfreqs(fasta):
	freqs = [] #[afreq, gfreq, cfreq, ufreq]
	a = 0
	u = 0
	c = 0
	g = 0
	tot = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seq = str(record.seq.transcribe().upper())
		a += seq.count('A')
		u += seq.count('U')
		c += seq.count('C')
		g += seq.count('G')
		tot += len(seq)

	freqs = [a/float(tot), g/float(tot), c/float(tot), u/float(tot)]

	return freqs

def getfreqs_boxplot(fasta, outfile, classid):
	freqs = {} # {seqname : [A,G,C,U]}
	for record in SeqIO.parse(fasta, 'fasta'):
		seq = str(record.seq.transcribe().upper())
		seqname = record.id
		tot = float(len(seq))
		if tot < 100:
			continue
		a = seq.count('A') / tot
		u = seq.count('U') / tot
		c = seq.count('C') / tot
		g = seq.count('G') / tot
		freqs[seqname] = [str(a),str(g),str(c),str(u)]

	outfh = open(outfile, 'w')
	outfh.write(('\t').join(['A','G','C','U','Class']) + '\n')
	for seq in freqs:
		outfh.write(('\t').join([freqs[seq][0], freqs[seq][1], freqs[seq][2], freqs[seq][3], classid]) + '\n')
	outfh.close()

#outfh = open(sys.argv[2], 'w')
#freqs = getfreqs(sys.argv[1])
#freqs = [str(freq) for freq in freqs]
#outfh.write(('\t').join(freqs) + '\t' + sys.argv[3] + '\n')
#outfh.close()

getfreqs_boxplot(sys.argv[1], sys.argv[2], sys.argv[3])