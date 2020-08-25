#Required modules: biopython, numpy, viennarna
#Given a fasta containing one sequence and a kmer you are interested in, make a fasta containing mutated sequences
#with the mutations occuring at a given rate everywhere except where the kmer occurs.  Then, using RNAfold, calculate
#the average basepair probabilites for bases in the original motif.

#Usage: python RNAfoldsimulation_motifbpprob.py -h


import sys
import subprocess
from Bio import SeqIO
import random
import argparse
import re


def getsummedbpprobs(seq, kmer, motifpositions):
	#Given a sequence, make bp probabilities and return the summed bp probability for a base (the sum of its pairing to every other base)
	probs = {} # {base : bpprob}

	#Populate dictionary
	for i in range(len(seq)):
		probs[i] = 0
	
	command = 'RNAfold -p'
	job = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	job.stdin.write(seq)
	output = job.communicate()
	structure = str(output[0].split('\n')[1].split(' ')[0])
	#print structure
	bpfh = open('dot.ps', 'r')
	#lines containing bp probs have 'ubox' in the 4th field
	for line in bpfh:
		line = line.strip().split(' ')
		if len(line) != 4:
			continue
		if line[3] == 'ubox':
			leftbase = int(line[0]) - 1 #make these 0 based
			rightbase = int(line[1]) - 1
			#3rd field is square root of bpprob
			bpprob = float(line[2])**2
			probs[leftbase] += bpprob
			probs[rightbase] += bpprob

	bpfh.close()
	
	#Get the bp prob over the motif and 3 nt up/down stream
	motifprobs = []
	motifprobs.append(probs[motifpositions[0] - 3])
	motifprobs.append(probs[motifpositions[0] - 2])
	motifprobs.append(probs[motifpositions[0] - 1])
	for i in range(len(kmer)):
		motifprobs.append(probs[motifpositions[i]])
	motifprobs.append(probs[motifpositions[-1] + 1])
	motifprobs.append(probs[motifpositions[-1] + 2])
	motifprobs.append(probs[motifpositions[-1] + 3])

	return motifprobs


def makemutatedfasta_constkmer(inputfasta, iterations, kmer):
	#Make a mutated sequence, but if you see "kmer", don't mutate it
	mutdict = {'A' : ['C','U','G'], 'G' : ['A','U','C'], 'C' : ['A','U','G'], 'U' : ['A','C','G'], 'T' : ['A','C','G']}
	lottery = []
	outfh = open('MutatedSequences.fasta', 'w')
	for i in range(1000):
		lottery.append(i + 1)

	#Get the original sequence and motif positions within it
	for record in SeqIO.parse(inputfasta, 'fasta'):
		seq = str(record.seq.transcribe())
		mp = [m.start() for m in re.finditer('(?={0})'.format(kmer), seq)] #list of all motif positions (motif start, 0 based)
		kmerlength = len(kmer)
		motifpositions = [] #list of all motif positions (any position covered by kmer)
		
		for motifpos in mp:
			for i in range(kmerlength):
				motifpositions.append(motifpos + i)

	for i in range(int(iterations)):
		if (i + 1) % 10000 == 0:
			print 'Creating mutated sequence {0} of {1}...'.format(i + 1, iterations)
		seqtitle = 'Seq' + str(i + 1)
		currentseq = ''
		for idx, nt in enumerate(seq):
			if idx in motifpositions: #if this position is part of a motif
				currentseq += nt
			else:
				ticket = random.choice(lottery)
				if ticket <= 940:
					currentseq += nt
				elif ticket >= 941 and ticket <= 960:
					mutation = mutdict[nt][0]
					currentseq += mutation
				elif ticket >= 961 and ticket <= 980:
					mutation = mutdict[nt][1]
					currentseq += mutation
				elif ticket >= 981 and ticket <= 1000:
					mutation = mutdict[nt][2]
					currentseq += mutation

		outfh.write('>' + seqtitle + '\n' + currentseq + '\n')

	print 'Done making mutated fasta!!'

	outfh.close()

	return motifpositions

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Original sequence in fasta format.')
	parser.add_argument('--kmer', type = str, help = 'RNA sequence of kmer of interest. Will not be mutated.')
	parser.add_argument('--iterations', type = int, help = 'Number of mutated sequences to make.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	outfh = open(args.outfile, 'w')
	motifpositions = makemutatedfasta_constkmer(args.fasta, args.iterations, args.kmer)
	
	#Extended motif positions are 3 bp around motif
	extendedmotifpositions = [motifpositions[0] - 3, motifpositions[0] - 2, motifpositions[0] - 1]
	for pos in motifpositions:
		extendedmotifpositions.append(pos)
	extendedmotifpositions.append(motifpositions[-1] + 1)
	extendedmotifpositions.append(motifpositions[-1] + 2)
	extendedmotifpositions.append(motifpositions[-1] + 3)

	outfh.write('Seqname' + '\t' + ('\t').join([str(pos) for pos in extendedmotifpositions]) + '\n')

	foldedseqcount = 0
	for record in SeqIO.parse('MutatedSequences.fasta', 'fasta'):
		foldedseqcount +=1
		if foldedseqcount % 10000 == 0:
			print 'Folding sequence {0}...'.format(foldedseqcount)
		seqname = record.id
		seq = str(record.seq)
		bpprobaroundmotif = getsummedbpprobs(seq, args.kmer, motifpositions)
		outfh.write(seqname + '\t' + ('\t').join([str(prob) for prob in bpprobaroundmotif]) + '\n')

	outfh.close()