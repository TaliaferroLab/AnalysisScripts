#python3

import subprocess
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import argparse
import os
import time

def fold_a_fasta_mfe(fasta):
	#Given a chopped fasta seqs, fold each chopped seq and count the number of gquads
	#Also, get the median folding energy over the whole fasta (which corresponds to one original oligo)
	#This fasta is designed to be a chopped fasta where every seq is a substring of a single oligo
	fragdict = {} #{fragid : [seqlength, energy, hasgquad]}
	gquadcounts = []
	energies = []
	command = 'RNAfold -g < {0}'.format(fasta)
	job = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE,
		stderr = subprocess.PIPE, stdin = subprocess.PIPE, encoding='utf8')
	output = job.communicate()[0][1:] #take first value of tuple and remove first character ('>')
	perseq = output.split('>')
	for seq in perseq:
		seqname = seq.split('\n')[0]
		ntseq = seq.split('\n')[1]
		structure = seq.split('\n')[2].split(' ')[0]
		energy = seq.split('\n')[2].split(' (')[1].replace(' ', '')[:-1]
		energy = float(energy)
		energies.append(energy)
		gquadcount = structure.count('+')
		gquadcounts.append(gquadcount)
		if gquadcount > 0:
			hasgquad = 'yes'
		else:
			hasgquad = 'no'
		fragdict[seqname] = [len(ntseq), energy, hasgquad]
	
	#Metrics that deal with the entire fasta (the entire oligo)
	medianenergy = np.median(energies)
	maxenergy = min(energies)
	if any(gquadcount > 0 for gquadcount in gquadcounts):
		hasgquad = 'yes'
	else:
		hasgquad = 'no'

	for f in os.listdir('.'):
		if f.endswith('.ps'):
			os.remove(f)

	return fragdict, medianenergy, maxenergy, hasgquad

#Parse the output of an RNA fold (usually the folding of one of several sequences)
#this is needed for mode == bpprob
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


def fold_a_fasta_bpprob(fasta):
	bpprobs = {} #{position (0-based) : [list of bpprobs for this position]}
	command = 'RNAfold -p -g < {0} > /dev/null'.format(fasta)
	subprocess.Popen(command, shell = True)

	#Now parse the RNAfold outputs
	for record in SeqIO.parse(fasta, 'fasta'):
		rnafoldoutputfile = (str(record.id) + '_dp.ps').replace('|', '_')
		windowstart = int(str(record.id).split('.')[2]) 
		probs = parseRNAfoldoutput(rnafoldoutputfile, str(record.seq))

		for position in probs:
			realseqposition = windowstart + position
			if realseqposition not in bpprobs:
				bpprobs[realseqposition] = [probs[position]]
			else:
				bpprobs[realseqposition].append(probs[position])

	#Take medians of the bp probs at every position across all windows as the bpprob for that position
	medianbpprobs = {}
	for position in bpprobs:
		medianbpprobs[position] = np.median(bpprobs[position])

	#Cleanup
	for f in os.listdir('.'):
		if f.endswith('_ss.ps') or f.endswith('_dp.ps'):
			os.remove(os.path.abspath(f))

	return medianbpprobs


def chopfasta(seqname, seq, windowsize, slidesize):
	#Given an seq, chop it into smaller seqs of length fragsize, sliding slidesize across the sequence
	#Make a fasta record out of it, then combine all of these fasta records into one file
	currentlocation = 0
	fastarecords = []
	#Remove adapters (first 20 nt and last 20 nt of sequence)
	seq = seq[20:-20]
	if len(seq) >= windowsize:
		while currentlocation + windowsize <= len(seq):
			currentseq = seq[currentlocation : currentlocation + windowsize]
			fastarecord = SeqRecord(Seq(currentseq, IUPAC.IUPACAmbiguousDNA), id = '{0}.{1}'.format(seqname, currentlocation), description = '')
			fastarecords.append(fastarecord)
			currentlocation += slidesize
		SeqIO.write(fastarecords, 'currentseqs.fasta', 'fasta')

def iteratefasta_mfe(fasta, fragoutput, oligooutput, windowsize, slidesize):
	#Get total number of seqs in fasta
	totalseqs = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		totalseqs +=1

	seqcounter = 0

	#Prepare outputs
	with open(fragoutput, 'w') as f:
		f.write(('\t').join(['fragname', 'oligoname', 'seqlength', 'energy', 'hasgquad']) + '\n')
	with open(oligooutput, 'w') as f:
		f.write(('\t').join(['oligoname', 'medianenergy', 'maxenergy', 'hasgquad']) + '\n')

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 100 == 0:
			print('Analyzing sequence {0} of {1}...'.format(seqcounter, totalseqs))
		if len(str(record.seq)) < windowsize:
			continue
		chopfasta(str(record.id), str(record.seq.transcribe()), windowsize, slidesize)
		fragdict, medianenergy, maxenergy, oligohasgquad = fold_a_fasta_mfe('currentseqs.fasta')
		with open(fragoutput, 'a') as f:
			for frag in fragdict:
				seqlength = fragdict[frag][0]
				energy = fragdict[frag][1]
				fraghasgquad = fragdict[frag][2]
				f.write(('\t').join([frag, ('.').join(frag.split('.')[:2]), str(seqlength), str(energy), str(fraghasgquad)]) + '\n')

		with open(oligooutput, 'a') as f:
			f.write(('\t').join([str(record.id), str(medianenergy), str(maxenergy), oligohasgquad]) + '\n')

		os.remove('currentseqs.fasta')

def iteratefasta_bpprob(fasta, oligooutput, windowsize, slidesize):
	#Get total number of seqs in fasta
	totalseqs = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		totalseqs +=1

	seqcounter = 0

	#Prepare outputs
	with open(oligooutput, 'w') as f:
		f.write(('\t').join(['oligoname', 'position', 'medianbpprob']) + '\n')

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 100 == 0:
			print('Analyzing sequence {0} of {1}...'.format(seqcounter, totalseqs))
		if len(str(record.seq)) < windowsize:
			continue
		chopfasta(str(record.id), str(record.seq.transcribe()), windowsize, slidesize)
		medianbpprobs = fold_a_fasta_bpprob('currentseqs.fasta')

		with open(oligooutput, 'a') as f:
			for pos in sorted(medianbpprobs):
				f.write(('\t').join([str(record.id), str(pos), str(medianbpprobs[pos])]) + '\n')

		os.remove('currentseqs.fasta')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to fold.')
	parser.add_argument('--fragoutput', type = str, help = 'Fragment output file.')
	parser.add_argument('--oligooutput', type = str, help = 'Oligo output file.')
	parser.add_argument('--windowsize', type = int, help = 'Size of oligo fragment to fold.')
	parser.add_argument('--slidesize', type = int, help = 'Distance to slide window along oligo.')
	parser.add_argument('--mode', type = str, choices = ['MFE', 'bpprob'], help = 'MFE or bpprob')
	args = parser.parse_args()

	if args.mode == 'MFE':
		iteratefasta_mfe(args.fasta, args.fragoutput, args.oligooutput, args.windowsize, args.slidesize)
	elif args.mode == 'bpprob':
		iteratefasta_bpprob(args.fasta, args.oligooutput, args.windowsize, args.slidesize)

