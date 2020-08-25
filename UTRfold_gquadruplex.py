#python2

import subprocess
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import argparse
import os


def fold_a_seq_gquad(seq):
	"""
	Returns number of G's predicted to be in quadruplexes by RNAfold
	Not really used here.
	"""

	command = 'RNAfold -g'
	job = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE,
		stderr = subprocess.PIPE, stdin = subprocess.PIPE)
	job.stdin.write(seq)
	output = job.communicate()
	structure, energy = (output[0].split('\n')[1].split(' ')[0], float(output[0].replace('(', ' ').replace(')', ' ').split()[-1])) #get secondary structure
	gquadcount = structure.count('+')
	return gquadcount, energy

def fold_a_fasta_gquad(fasta):
	#Given a chopped fasta seqs, fold each chopped seq and count the number of gquads
	#Also, get the median folding energy over the whole fasta (which corresponds to one original seq)
	gquads = 0
	gquadpositions = []
	energies = []
	command = 'RNAfold -g < {0}'.format(fasta)
	job = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE,
		stderr = subprocess.PIPE, stdin = subprocess.PIPE)
	output = job.communicate()[0][1:] #take first value of tuple and remove first character ('>')
	perseq = output.split('>')
	for seq in perseq:
		startcoord = int(seq.split('\n')[0].split('__')[-1]) #0 based start coord of this fragment of the fasta
		structure = seq.split('\n')[2].split(' ')[0]
		energy = seq.split('\n')[2].split(' (')[1].replace(' ', '')[:-1]
		gquads += structure.count('+')
		#Get positions of quadruplexed G residues with respect to their coordinate in the original (unchopped) sequence. 0-based.
		gquadpos = [pos + startcoord for pos, char in enumerate(structure) if char == '+']
		energy = float(energy)
		energies.append(energy)
		gquadpositions += gquadpos
	
	medianenergy = np.median(energies)
	#Since this is a sliding window that seems some positions twice, there may be duplicates in gquadpositions
	gquadpositions = list(set(gquadpositions))
	#gquads may also be duplicated due to sliding window
	gquads = len(gquadpositions)

	for f in os.listdir('.'):
		if f.endswith('.ps'):
			os.remove(f)

	return gquads, medianenergy, gquadpositions


def chopfasta(seqname, seq, windowsize, slidesize):
	#Given a seq, chop it into smaller fasta seqs of length fragsize, sliding slidesize across the sequence
	#Make a fasta record out of it, then combine all of these fasta records into one file
	currentlocation = 0
	fastarecords = []
	if len(seq) >= windowsize:
		while currentlocation + windowsize <= len(seq):
			currentseq = seq[currentlocation : currentlocation + windowsize]
			fastarecord = SeqRecord(Seq(currentseq, IUPAC.IUPACAmbiguousDNA), id = 'seqname__{0}'.format(currentlocation), description = '')
			fastarecords.append(fastarecord)
			currentlocation += slidesize
		SeqIO.write(fastarecords, 'currentseqs.fasta', 'fasta')


def iteratefasta(fasta, output, windowsize, slidesize):
	#Get total number of seqs in fasta
	totalseqs = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		totalseqs +=1

	seqcounter = 0

	#Prepare output
	with open(output, 'w') as f:
		f.write(('\t').join(['seqname','seqlength','gquadcount','medianenergy','gquadpos']) + '\n')

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 100 == 0:
			print 'Analyzing sequence {0} of {1}...'.format(seqcounter, totalseqs)
		if len(str(record.seq)) < windowsize:
			continue
		chopfasta(str(record.id), str(record.seq.transcribe()), windowsize, slidesize)
		gquads, medianenergy, gquadpositions = fold_a_fasta_gquad('currentseqs.fasta')
		if gquadpositions:
			gquadpositions = ','.join([str(gquadpos) for gquadpos in sorted(gquadpositions)])
		else:
			gquadpositions = 'none'
		os.remove('currentseqs.fasta')

		with open(output, 'a') as f:
			f.write(('\t').join([str(record.id), str(len(str(record.seq))), str(gquads), str(medianenergy), gquadpositions]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to fold.')
	parser.add_argument('--output', type = str, help = 'Output file.')
	parser.add_argument('--windowsize', type = int, help = 'Size of UTR fragment to fold.')
	parser.add_argument('--slidesize', type = int, help = 'Distance to slide window along UTR.')
	args = parser.parse_args()
	iteratefasta(args.fasta, args.output, args.windowsize, args.slidesize)


