#Starting with a gff of regions you are interested in, define new regions that are windows +/- 40 bp of every
#nt in the gff.  Each window is associated with a kmer which is positions 41 to 41 + k within 
#that window, i.e. the kmer starts at the nucleotide you are at during your walk through every
#nt in the gff.  Then fold that window using RNAstructure.  Record the average bp probability over the kmer.
#Make a bed of these kmer positions (which is every possible kmer in your original gff) and pass that to Jason's
#branch length script.  Then try to correlate the branch length of a particular kmer with its bp probability.

#Output is in kmerbpprobs.bed

from Bio import SeqIO
import sys
import subprocess
import os
from numpy import mean
import argparse

def getregionstofold(gff, k):
	#Starting with a gff (of oligos, for example), make gff lines where each line is +/- 40 bp from beginning of every kmer
	#in the original gff line
	k = int(k)
	print 'Identifying regions to fold...'
	rtfgff = [] #"Regions to fold gff"
	gfffh = open(gff, 'r')
	for line in gfffh:
		line = line.strip().split('\t')
		chrm = line[0]
		species = line[1]
		location = line[2]
		start = int(line[3])
		stop = int(line[4])
		strand = line[6]
		ID = line[8][3:]

		#Want to exclude the first 10 nt of the 5' portion of introns and the last 20 nt of the 3' portion of introns
		#to exclude anything near splice sites.
		if strand == '+':
			if location == 'up_5intron' or location == 'dn_5intron':
				windowstart = 10
				windowstop = 0
			elif location == 'up_3intron' or location == 'dn_3intron':
				windowstart = 0
				windowstop = 20
		elif strand == '-':
			if location == 'up_5intron' or location == 'dn_5intron':
				windowstart = 0
				windowstop = 10
			elif location == 'up_3intron' or location == 'dn_3intron':
				windowstart = 20
				windowstop = 0

		for i in range(start + windowstart, stop - k - windowstop + 1):
			kmerstart = i
			kmerstop = i + k
			rtfstart = i - 40
			rtfstop = i + 40
			#new ID is original ID with kmer coords attached after '_'
			newID = 'ID=' + ID + ';' + str(kmerstart) + '|' + str(kmerstop)
			rtfgff.append([chrm, species, location, rtfstart, rtfstop, '.', strand, '.', newID])

	gfffh.close()
	return rtfgff

#The output from ProbabilityPlot is the probability of all possible basepairs. I just want to know
#what's the probability of a base to be paired to any other base.  So, for each nucleotide, sum all of
#the basepair probabilities for that nucleotide.
def summarizebpprobmatrix(bpprobmatrix):
    infh = open(bpprobmatrix, 'r')
    bpprobs = {} # {nt (a number, its 1-based position in the folded sequence) : cumulative_bpprob}
    for line in infh:
        line = line.strip().split('\t')
        if len(line) != 3 or line[0] == 'i': #skip header lines
            continue
        nt1 = int(line[0])
        nt2 = int(line[1])
        prob = 10**(float(line[2]) * -1)
        if nt1 not in bpprobs:
            bpprobs[nt1] = prob
        elif nt1 in bpprobs:
            newprob = bpprobs[nt1] + prob
            bpprobs[nt1] = newprob
        if nt2 not in bpprobs:
            bpprobs[nt2] = prob
        elif nt2 in bpprobs:
            newprob = bpprobs[nt2] + prob
            bpprobs[nt2] = newprob

    infh.close()

    return bpprobs

def foldrtf(rtfgff, genomefasta, k):
	#rtf = region to fold
	k = int(k)
	print 'Indexing genome...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	print 'Done indexing!'
	counter = 0
	bedfh = open('kmerbpprobs.bed', 'w')

	for region in rtfgff:
		counter +=1
		if counter % 100 == 0:
			print '\n***\nFolding region {0} of {1}...\n***\n'.format(counter, len(rtfgff))
		
		chrm, species, location, rtfstart, rtfstop, strand, ID = region[0], region[1], region[2], region[3], region[4], region[6], region[8]
		kmerstart = ID.split(';')[1].split('|')[0]
		kmerstop = ID.split(';')[1].split('|')[1]
		if strand == '+':
			seq = seq_dict[chrm].seq[rtfstart:rtfstop].transcribe().upper()
		elif strand == '-':
			seq = seq_dict[chrm].seq[rtfstart:rtfstop].reverse_complement().transcribe().upper()

		outfh = open('temp.fasta', 'w')
		outfh.write('>' + ID + '\n' + str(seq))
		outfh.close()

		#Make the partition file
		subprocess.call(['partition', 'temp.fasta', 'output.pfs'])
		#Make the bp prob matrix
		subprocess.call(['ProbabilityPlot', '--text', 'output.pfs', 'bpprobmatrix.txt'])

		bpprobs = summarizebpprobmatrix('bpprobmatrix.txt')
		kmerbpprobs = []
		
		#If we are on the + strand, the kmer of interest begins at position 41 (1 based, which is what the bpprob file is)
		if strand == '+':
			for i in range(41,41 + k):
				try:
					kmerbpprobs.append(bpprobs[i])
				except KeyError:
					kmerbpprobs.append(0)
		#If we are on the - strand, the kmer of interest begins at position 41-k (1 based)
		elif strand == '-':
			for i in range(41-k, 41):
				try:
					kmerbpprobs.append(bpprobs[i])
				#If there were any nt with 0 bp probability, they dont get recorded
				except KeyError:
					kmerbpprobs.append(0)

		meankmerbpprob = mean(kmerbpprobs)
		
		#Make a new bed (not gff) line where the mean kmer bpprob is in the 6th field (index 5) and the coordinates are for the kmer
		rtfbed_bpprob = [chrm, str(kmerstart), str(kmerstop), ID, str(meankmerbpprob), strand]
		print ('\t').join(rtfbed_bpprob)
		with open('kmerbpprobs.bed', 'a') as bedfile:
			bedfile.write(('\t').join(rtfbed_bpprob) + '\n')
		os.remove('output.pfs')
		os.remove('bpprobmatrix.txt')



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF file of introns to consider.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('-k', type = int, help = 'Size of kmers.')
	args = parser.parse_args()

	rtfgff = getregionstofold(args.gff, args.k)
	foldrtf(rtfgff, args.genomefasta, args.k)