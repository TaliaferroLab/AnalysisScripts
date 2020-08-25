#Take a fasta of sequences.  Fold it using a sliding window approach and return base pair probabilities for every nt.
#Each nt may or may not have multiple values assigned to it, depending on how many windows it fell in to.
#Take the mean bp prob across all windows as the bp prob for that nt.
#Then, divide into bins and take the mean bp of all nt in that bin.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import subprocess
import time
import os
import sys
import argparse



def chopfasta(seq, windowsize, slidesize):
	#Given a seq, chop it into smaller fasta seqs of length fragsize, sliding slidesize across the sequence
	#Make a fasta record out of it, then combine all of these fasta records into one file
	currentlocation = 0
	fastarecords = []
	if len(seq) >= windowsize:
		while currentlocation + windowsize <= len(seq):
			currentseq = seq[currentlocation : currentlocation + windowsize]
			fastarecord = SeqRecord(Seq(currentseq, IUPAC.IUPACAmbiguousDNA), id = 'seq__{0}'.format(currentlocation), description = '')
			fastarecords.append(fastarecord)
			currentlocation += slidesize

	#Now start at the end and go backwards to ensure that every nt is covered in at least one window
		currentlocation = len(seq) - windowsize
		while currentlocation >= 0:
			currentseq = seq[currentlocation : currentlocation + windowsize]
			fastarecord = SeqRecord(Seq(currentseq, IUPAC.IUPACAmbiguousDNA), id = 'seq__{0}'.format(currentlocation), description = '')
			fastarecords.append(fastarecord)
			currentlocation -= slidesize

		SeqIO.write(fastarecords, 'currentseqs.fasta', 'fasta')



#Parse the output of an RNA fold (usually the folding of one of several sequences)
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

#Parse the output of an RNA fold (usually the folding of one of several sequences)
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

#Given a sequence (of one CDS or one UTR, for example), chop it into windowed sizes.
#Then fold each window.
#Keep track of the bp prob for each nt in each window.

def foldseq(seq, windowsize, slidesize):
	#If the sequence is smaller than the windowsize, we can't deal with it
	if len(seq) < windowsize:
		return None

	bpprobs = {} #{position (0-based) : [list of bpprobs for this position]}
	
	#populate bpprobs
	for i in range(len(seq)):
		bpprobs[i] = []

	#Chop the sequence into windows
	chopfasta(seq, windowsize, slidesize)

	#Fold the fasta file, send stdout to trash
	command = 'RNAfold -p < currentseqs.fasta > /dev/null'
	subprocess.Popen(command, shell = True)

	#Now parse the RNAfold outputs
	for record in SeqIO.parse('currentseqs.fasta', 'fasta'):
		rnafoldoutputfile = str(record.id) + '_dp.ps'
		windowstart = int(str(record.id).split('__')[1])
		probs = parseRNAfoldoutput(rnafoldoutputfile, str(record.seq))

		for position in probs:
			realseqposition = windowstart + position
			bpprobs[realseqposition].append(probs[position])

	#Take medians of the bp probs at every position across all windows as the bpprob for that position
	medianbpprobs = {}
	for position in bpprobs:
		medianbpprobs[position] = np.median(bpprobs[position])

	#Cleanup
	for f in os.listdir('.'):
		if f.endswith('_ss.ps') or f.endswith('_dp.ps'):
			os.remove(os.path.abspath(f))
	os.remove('currentseqs.fasta')

	return medianbpprobs


#Now turn the medianbpprobs into bins (like a metagene)
def turnintobins(medianbpprobs, numberofbins):
	#Given a list of nt positions, split it into n roughly-equal parts
	#http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
	listofpositions = medianbpprobs.keys()
	bins = []
	k, m = divmod(len(listofpositions), numberofbins)
	gen = (listofpositions[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(numberofbins))
	for b in gen:
		bins.append(b)

	binnedbpprobs = {} #{bin : [bp probs of nt in that bin]}

	for ind, b in enumerate(bins):
		binnumber = ind + 1
		binnedbpprobs[binnumber] = []
		for nt in b:
			binnedbpprobs[binnumber].append(medianbpprobs[nt])

	#Now take the mean of all bpprobs in each bin as the bpprob for that bin
	meanbinnedbpprobs = {} #{bin : mean bpprob of nt in that bin}
	for binnumber in binnedbpprobs:
		meanbinnedbpprobs[binnumber] = np.mean(binnedbpprobs[binnumber])

	return meanbinnedbpprobs


#Now iterate over an entire fasta of sequences
def iteratefasta(fasta, windowsize, slidesize, numberofbins, outputfile, region, seqclass):
	if os.path.isfile(outputfile) == False:
		with open(outputfile, 'w') as outfh:
			outfh.write(('\t').join(['seqname', 'binnumber', 'bpprob', 'region', 'seqclass']) + '\n')
	
	for record in SeqIO.parse(fasta, 'fasta'):
		seq = str(record.seq.transcribe().upper())
		seqname = str(record.id)
		medianbpprobs = foldseq(seq, windowsize, slidesize)
		
		if not medianbpprobs: #if the length of the sequence is less than windowsize, medianbpprobs returns None
			continue
			
		meanbinnedbpprobs = turnintobins(medianbpprobs, numberofbins)
		with open(outputfile, 'a') as outfh:
			for binnumber in meanbinnedbpprobs:
				meanbpprob = meanbinnedbpprobs[binnumber]
				outfh.write(('\t').join([seqname, str(binnumber), str(meanbpprob), region, seqclass]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to fold.')
	parser.add_argument('--windowsize', type = int, help = 'Width of subsequences that get folded.')
	parser.add_argument('--slidesize', type = int, help = 'How far the window gets slid along the sequence to define the next window.')
	parser.add_argument('--numberofbins', type = int, help = 'How many bins to use for the metagene.')
	parser.add_argument('--region', type = str, help = 'Transcript region these sequences came from.  UTR? CDS?')
	parser.add_argument('--seqclass', type = str, help = 'Class of these sequences.  RO up? LR down?')
	parser.add_argument('--outputfile', type = str, help = 'Output file.')

	args = parser.parse_args()
	iteratefasta(args.fasta, args.windowsize, args.slidesize, args.numberofbins, args.outputfile, args.region, args.seqclass)
