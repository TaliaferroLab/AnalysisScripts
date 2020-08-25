#Given a fasta (probably of CDS sequences), determine the fraction of instances of a particular 3mer
#that are in the +1 (coding) frame.

#Usage: python kmersinframe fasta kmer label


from Bio import SeqIO
import sys

def getframefractions(fasta, kmer_of_interest):
	#Count number of kmer_of_interest in each sequence.
	#See what fraction of all kmers are that kmer of interest.
	k = len(kmer_of_interest)
	kmers = {} #{seqname : [kmersofinterest in frame 1, kmersofinterest in frame 2, kmersofinterest in frame 3]}
	frame1fracs = [] #fractions (1 per seq) of kmers of interest that are in frame 1

	for record in SeqIO.parse(fasta, 'fasta'):
		frame1 = 0
		frame2 = 0
		frame3 = 0
		totalkmersofinterest = 0

		seq = str(record.seq.transcribe())

		for i in range(len(seq) - k + 1):
			kmer = seq[i:i+k]
			frame = (i + 1) % 3
			if frame == 0:
				frame = 3
			if kmer == kmer_of_interest:
				if frame == 1:
					frame1 +=1
				elif frame == 2:
					frame2 +=1
				elif frame == 3:
					frame3 +=1

				totalkmersofinterest +=1

		kmers[record.id] = [frame1, frame2, frame3]
		if totalkmersofinterest > 0:
			frame1frac = frame1 / float(totalkmersofinterest)
			frame1fracs.append(frame1frac)

	
	return frame1fracs

def getleucinerepeatfreq(fasta):
	#Look for how often leucine rich repeats occur in a fasta of CDS sequence (nt seq)
	#Leucine rich repeat seq = LxxLxLxxNxL L= L,I,V,P, N = N,T,C,S
	k = 11 #length of LRR
	L = ['L','I','V','P']
	N = ['N','T','C','S']
	totalleucinecount = 0
	totalLRRcount = 0
	totalAAcount = 0

	for record in SeqIO.parse(fasta, 'fasta'):
		AAseq = str(record.seq.translate())
		leucinecount = AAseq.count('L')
		totalleucinecount += leucinecount
		totalAAcount += len(AAseq)
		LRRcount = 0

		for i in range(len(AAseq) - k + 1):
			kmer = AAseq[i:i+k]
			if kmer[0] in L and kmer[3] in L and kmer[5] in L and kmer[8] in N and kmer[10] in L:
				LRRcount += 1

		totalLRRcount += LRRcount
	
	print totalLRRcount / float(totalleucinecount)

		


frame1fracs = getframefractions(sys.argv[1], sys.argv[2])

for frac in frame1fracs:
	print str(frac) + '\t' + sys.argv[3]

#getleucinerepeatfreq(sys.argv[1])