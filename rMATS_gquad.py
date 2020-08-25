from Bio import SeqIO
import subprocess
import os
import pandas as pd
import gzip
import sys


#Given a sequence, return the positions that are participating in a quadruplex in structures
#that are within 2 kcal/mol of the MFE.
def getgquadpos(seq):
	gquadpos = [False] * len(seq)
	ps = subprocess.Popen(['echo', seq], stdout = subprocess.PIPE)
	output = subprocess.check_output(['RNAsubopt', '-g', '-e', '2'], stdin = ps.stdout)
	#parse output
	output = output.split('\n')[1:]
	for fold in output:
		struct = fold.split(' ')[0]
		for pos, char in enumerate(struct):
			if char == '+':
				gquadpos[pos] = True

	gquadpos = [pos for pos, char in enumerate(gquadpos) if char == True]
	return gquadpos

def getSequences(rMATStable, genomefasta):
	#Make (or blank) files
	for eventclass in ['moreincluded', 'lessincluded', 'control']:
		for seqregion in ['seq1', 'seq2', 'seq3', 'seq4']:
			outfh = open('{0}.{1}.fasta'.format(eventclass, seqregion), 'w')
			outfh.close()

	#Index fasta
	print 'Indexing genome...'
	if os.path.basename(genomefasta).endswith('.gz'):
		seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	else:
		seq_dict = SeqIO.to_dict(SeqIO.parse(open(genomefasta, 'r'), 'fasta'))
	print 'Done indexing!'

	df = pd.read_table(rMATStable)

	#|||||||||---------------//-------------||||||||//|||||||||-----------------//-----------------||||||||||||||||| (this is on the plus strand)
	# seq1 -50 to + 150           seq2 -150 to + 50         seq3 -50 to + 150             seq4 -150 to + 50
	#            <stretch1>     <stretch2>                              <stretch3>          <stretch4>

	for index, row in df.iterrows():
		chrm = row['chr']
		strand = row['strand']
		exonstart = row['exonStart_0base']
		exonend = row['exonEnd']
		upstreamexonstart = row['upstreamES']
		upstreamexonend = row['upstreamEE']
		downstreamexonstart = row['downstreamES']
		downstreamexonend = row['downstreamEE']
		eventID = chrm + ':' + strand + ':' + str(upstreamexonstart) + '-' + str(upstreamexonend) + ':' + str(exonstart) + '-' + str(exonend) + ':' + str(downstreamexonstart) + '-' + str(downstreamexonend)
		FDR = row['FDR']
		dpsi = row['IncLevelDifference']

		if row['FDR'] < 0.05 and row['IncLevelDifference'] > 0:
			eventclass = 'moreincluded'
		elif row['FDR'] < 0.05 and row['IncLevelDifference'] < 0:
			eventclass = 'lessincluded'
		elif row['FDR'] > 0.2:
			eventclass = 'control'
		else:
			continue

		#Control for +150 or -150 sequences stretching into the next exon
		#If it does, cut it off at the exon/intron boundary
		if strand == '+':
			stretch1 = min(150, exonstart - upstreamexonend)
			stretch2 = min(150, exonstart - upstreamexonend)
			stretch3 = min(150, downstreamexonstart - exonend)
			stretch4 = min(150, downstreamexonstart - exonend)

		elif strand == '-':
			stretch1 = min(150, downstreamexonstart - exonend)
			stretch2 = min(150, downstreamexonstart - exonend)
			stretch3 = min(150, exonstart - upstreamexonend)
			stretch4 = min(150, exonstart - upstreamexonend)


		if strand == '+':
			seq1 = seq_dict[chrm].seq[upstreamexonend - 50 : upstreamexonend + stretch1]
			seq2 = seq_dict[chrm].seq[exonstart - stretch2 : exonstart + 50]
			seq3 = seq_dict[chrm].seq[exonend - 50 : exonend + stretch3]
			seq4 = seq_dict[chrm].seq[downstreamexonstart - stretch4 : downstreamexonstart + 50]

		elif strand == '-':
			seq1 = seq_dict[chrm].seq[downstreamexonstart - stretch1 : downstreamexonstart + 50].reverse_complement()
			seq2 = seq_dict[chrm].seq[exonend - 50 : exonend + stretch2].reverse_complement()
			seq3 = seq_dict[chrm].seq[exonstart - stretch3 : exonstart + 50].reverse_complement()
			seq4 = seq_dict[chrm].seq[upstreamexonend - 50 : upstreamexonend + stretch4].reverse_complement()

		for fn, seqregion in zip(['seq1', 'seq2', 'seq3', 'seq4'], [str(seq1).upper(), str(seq2).upper(), str(seq3).upper(), str(seq4).upper()]):
			with open('{0}.{1}.fasta'.format(eventclass, fn), 'a') as outfh:
				outfh.write('>' + eventID + '\n' + seqregion + '\n')

def foldfasta(fasta, windowsize, slidesize):
	seqsinfile = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqsinfile +=1

	seqcounter = 0
	fastagquaddict = {} #{seqname : {0-basedposition : <0 or 1>}} 1 if was quadruplexed in some folded structure in some window.  Will not doublecount.
	countergquaddict = {} #{0-basedposition : [number of seqs with no gquad here, number of seqs with a gquad here]} This is an overall summary dict.
	
	#Populate countergquaddict with 0s
	for i in range(200):
		countergquaddict[i] = [0, 0]

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 20 == 0:
			print 'Sequence {0} of {1}...'.format(seqcounter, seqsinfile)
		gquaddict = {} #{0-basedposition : <0 or 1>]} 1 if was quadruplexed in some folded structure in some window.  Will not doublecount.
		seq = str(record.seq)

		#Populate gquaddict with 0s
		for i in range(len(seq)):
			gquaddict[i] = 0

		windowstart = 0
		while windowstart + windowsize <= len(seq):
			windowseq = seq[windowstart : windowstart + windowsize]
			gquadpos = getgquadpos(windowseq)
			#Adjust gquadpos for startposition of the window
			gquadpos = [pos + windowstart for pos in gquadpos]
			for pos in gquadpos:
				gquaddict[pos] = 1
			windowstart += slidesize

		for pos in gquaddict:
			if gquaddict[pos] == 1:
				countergquaddict[pos][1] +=1
			elif gquaddict[pos] == 0:
				countergquaddict[pos][0] +=1
			else:
				print 'ERROR'
		fastagquaddict[record.id] = gquaddict

	#Write gquadpos per seq out to file
	outfile = os.path.splitext(fasta)[0] + '.gquadpos.txt'
	with open(outfile, 'w') as outfh:
		for seqid in fastagquaddict:
			gquadpos = []
			for pos in fastagquaddict[seqid]:
				if fastagquaddict[seqid][pos] == 1:
					gquadpos.append(pos)
			gquadpos = sorted(gquadpos)
			gquadpos = [str(pos) for pos in gquadpos]
			gquadpos = (',').join(gquadpos)
			if not gquadpos:
				gquadpos = 'None'
			outfh.write(seqid + '\t' + gquadpos + '\n')

	print countergquaddict
	return countergquaddict

def iteratefastas():
	outfh = open('Gquadresults.txt', 'w')
	outfh.close()

	fastas = []
	for eventtype in ['lessincluded', 'control', 'moreincluded']:
		for region in ['seq1', 'seq2', 'seq3', 'seq4']:
			fastas.append('{0}.{1}.fasta'.format(eventtype, region))

	print fastas

	for fasta in fastas:
		countergquaddict = foldfasta(fasta, 60, 1)
		with open('Gquadresults.txt', 'a') as outfh:
			outfh.write(fasta)
			for pos in range(200):
				outfh.write('\t' + str(pos) + ',' + str(countergquaddict[pos][0]) + ',' + str(countergquaddict[pos][1])) #fasta pos,noquad,gquad pos,noquad,gquad
			outfh.write('\n')

#getSequences(sys.argv[1], sys.argv[2])
#foldfasta(sys.argv[1], 60, 1)

iteratefastas()

