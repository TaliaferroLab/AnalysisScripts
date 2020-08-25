import argparse
from Bio import SeqIO
import subprocess
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import sys
import os
import time

#Make all possible single mutations in a sequence
def makesinglemutations(seq, untouchablepositions, output):
	mutdict = {'A' : ['C','U','G'], 'G' : ['A','U','C'], 'C' : ['A','U','G'], 'U' : ['A','C','G'], 'T' : ['A','C','G']}
	seq = seq.replace('T', 'U')

	for i in range(len(seq)):
		if i in untouchablepositions: #don't make in mutations in the motif or adapters
			continue
		correctnt = seq[i]
		for mutation in mutdict[correctnt]:
			seqlist = list(seq)
			seqlist[i] = mutation
			mutatedseq = ''.join(seqlist)
			mutatedseqname = str(i + 1) + mutation
			with open(output, 'a') as f:
				f.write('>' + mutatedseqname + '\n' + mutatedseq + '\n')

def makepairwisemutations(seq, untouchablepositions, output):
	mutdict = {'A' : ['C','U','G'], 'G' : ['A','U','C'], 'C' : ['A','U','G'], 'U' : ['A','C','G'], 'T' : ['A','C','G']}
	seq = seq.replace('T', 'U')

	for i in range(len(seq)):
		if i in untouchablepositions:
			continue
		correcti = seq[i]
		for mutationi in mutdict[correcti]:
			seqlist = list(seq)
			seqlist[i] = mutationi
			for j in range(len(seq)):
				if j in untouchablepositions or i == j:
					continue
				correctj = seq[j]
				for mutationj in mutdict[correctj]:
					seqlist[j] = mutationj
					mutatedseq = ''.join(seqlist)
					mutatedseqname = str(i + 1) + mutationi + '_' + str(j + 1) + mutationj
					with open(output, 'a') as f:
						f.write('>' + mutatedseqname + '\n' + mutatedseq + '\n')

def parseRNAfoldoutput(outputfile, seq):
	#Given a sequence, make bp probabilities and return the summed bp probability for a base (the sum of its pairing to every other base)
	probs = {} # {base : bpprob}

	#Populate dictionary
	for i in range(len(seq)):
		probs[i] = 0

	#For some reason the slow point in the folding analysis is the writing of the .ps files.
	#This means that sometimes the folding will be finished, but the .ps file that this is trying to read is not made yet
	#Try to read it.  If it's not there, wait 1 sec and then try again.
	while True:
		try:
			bpfh = open(outputfile, 'r')
			break
		except IOError:
			time.sleep(1)

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

	return probs

def foldfasta(fasta, motifpositions, output):
	#Fold sequences in a set of muatated sequences in a fasta
	#Motif positions is a list of 1-based positions where the motif of interest is in these seqs.
	
	#See how many sequences there are in the fasta
	totalseqs = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		totalseqs +=1
		if totalseqs == 1:
			dummyseq = str(record.seq) #Use a sequence to feed to parseRNAfoldoutput so that it can get the lengths

	#Prepare output file
	with open(output, 'w') as f:
		f.write(('\t').join(['seqname', 'Pos1prob', 'Pos2prob', 'Pos3prob', 'Pos4prob', 'meanmotifbpprob']) + '\n')

	#Put 10000 seqs into a new fasta file
	fastarecords = []
	ids = []
	readcounter = 0
	analyzedseqs = 0 #total number of analyzed sequences
	for record in SeqIO.parse(fasta, 'fasta'):
		readcounter +=1
		analyzedseqs +=1
		ids.append(str(record.id))
		seq = str(record.seq)
		fastarecord = SeqRecord(Seq(seq, IUPAC.IUPACAmbiguousRNA), id = str(record.id), description = '')
		fastarecords.append(fastarecord)
		if readcounter == 10000 or analyzedseqs == totalseqs:
			#Write records to a temporary file
			SeqIO.write(fastarecords, 'currentreads.fasta', 'fasta')

			#Fold the fasta file, send stdout to trash
			command = 'RNAfold -p < currentreads.fasta > /dev/null'
			subprocess.Popen(command, shell=True)

			#Parse the outputs
			for ID in ids:
				probs = parseRNAfoldoutput('{0}_dp.ps'.format(ID), dummyseq)
				motifbpprobs = []
				for position in motifpositions:
					prob = probs[position - 1] #Positions are 0-based, motifpositions are 1-based
					motifbpprobs.append(prob)
				meanmotifbpprob = np.mean(motifbpprobs)
				with open(output, 'a') as f:
					f.write(('\t').join([ID, str(motifbpprobs[0]), str(motifbpprobs[1]), str(motifbpprobs[2]), str(motifbpprobs[3]), str(meanmotifbpprob)]) + '\n')

			#Reset stuff
			readcounter = 0
			fastarecords = []
			ids = []
			for f in os.listdir('.'):
				if f.endswith('_ss.ps') or f.endswith('_dp.ps'):
					os.remove(os.path.abspath(f))
			os.remove('currentreads.fasta')





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--sequence', type = str, help = 'Starting sequence.')
	parser.add_argument('--mutations', type = int, help = 'Number of mutations per molecule.')
	parser.add_argument('--untouchablepos', type = str, help = 'Comma-separated list of positions not to mutate (1-based, usually adapters and motifs), e.g. 1-23,133-150')
	parser.add_argument('--mode', type = str, choices = ['makemutations', 'foldseqs'], help = 'Do you want to make mutations or fold sequences?')
	parser.add_argument('--motifpositions', type = str, help = 'Comma-separated list of motif positions. 1-based. e.g. 56,57,58,59. Use only if in foldseqs mode.')
	parser.add_argument('--fasta', type = str, help = 'Fasta file of mutated seqs to fold.  Use only if in foldseqs mode.')
	parser.add_argument('--output', type = str, help = 'Output file.')
	args = parser.parse_args()

	if args.mode == 'makemutations':
		if args.motifpositions or args.fasta:
			print 'Cannot use motifpositions option or fasta option with makemutations mode.'
			sys.exit()
		untouchablepos = args.untouchablepos.split(',')
		untouchablepositions = []
		for posrange in untouchablepos:
			boundaries = posrange.split('-')
			for i in range(int(boundaries[0]), int(boundaries[1]) + 1):
				untouchablepositions.append(i - 1)

		if args.mutations == 1:
			makesinglemutations(args.sequence, untouchablepositions, args.output)

		elif args.mutations == 2:
			makepairwisemutations(args.sequence, untouchablepositions, args.output)

	elif args.mode == 'foldseqs':
		if args.sequence or args.mutations or args.untouchablepos:
			print 'Invalid options for use with mode foldseqs.  Only motifpositions, fasta, and output are necessary.'
			sys.exit()
		motifpositions = args.motifpositions.split(',')
		motifpositions = [int(position) for position in motifpositions]
		foldfasta(args.fasta, motifpositions, args.output)

