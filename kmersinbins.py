from Bio import SeqIO
import random
import argparse

def countKmer(fasta, kmer_of_interest):
	#Count number of kmer_of_interest in each sequence.
	#See what fraction of all kmers are that kmer of interest.
	k = len(kmer_of_interest)
	kmers = {} #{seqname : [kmer_of_interest count, total kmer count]}

	for record in SeqIO.parse(fasta, 'fasta'):
		kmercounts = 0
		totalkmercounts = 0

		seq = str(record.seq.transcribe())

		for i in range(len(seq) - k + 1):
			kmer = seq[i:i+k]
			totalkmercounts +=1
			if kmer == kmer_of_interest:
				kmercounts +=1

		kmers[record.id] = [kmercounts, totalkmercounts]

	return kmers

def splitfasta(fasta):
	#Split the fasta into 20 ~equally sized bins. Bin contents are randomly assigned so not every seq may be assigned a bin and some seqs may be assigned twice.
	seqnames = []
	seqbins = [] #list of seqnames in 20 bins (each bin itself is a list of names)
	for record in SeqIO.parse(fasta, 'fasta'):
		seqnames.append(record.id)

	no_of_seqs = len(seqnames)
	samplesize = no_of_seqs / 20
	for i in range(20):
		selected_seqs = random.sample(seqnames, 20)
		seqbins.append(selected_seqs)

	return seqbins

def analyzefasta(fasta, kmer_of_interest):
	#Intersect split fasta and kmer count dict
	binfreqs = {} # {binnumber : fraction of kmers that are kmer of interest}
	kmerdict = countKmer(fasta, kmer_of_interest)
	seqbins = splitfasta(fasta)
	i = 1
	for seqbin in seqbins:
		binname = 'bin{0}'.format(i)
		kmercount = 0
		totalkmercounts = 0
		for seq in kmerdict:
			if seq in seqbin:
				kmer_of_interest_count = kmerdict[seq][0]
				total_kmer_count = kmerdict[seq][1]
				kmercount += kmer_of_interest_count
				totalkmercounts += total_kmer_count

		kmerfraction = kmercount / float(totalkmercounts)
		binfreqs[binname] = kmerfraction
		i +=1

	return binfreqs



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--sequences', type = str, help = 'Comma separated list of fasta files to consider in the order that you want the fastas to be considered.')
	parser.add_argument('--kmersofinterest', type = str)
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	outfh = open(args.outfile, 'w')
	fastas = args.sequences.split(',')
	kmers = args.kmersofinterest.split(',')
	outfh.write(('\t').join(['kmer', 'Bin','subbin','kmerfreq']) + '\n')

	for kmer in kmers:
		fastacounter = 1
		for fasta in fastas:
			binname = 'Bin{0}'.format(fastacounter)
			binfreqs = analyzefasta(fasta, kmer)
			for b in sorted(binfreqs):
				outfh.write(('\t').join([kmer, binname, b, str(binfreqs[b])]) + '\n')
			fastacounter +=1

	outfh.close()




