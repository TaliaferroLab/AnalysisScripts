from Bio import SeqIO
import os
import pandas as pd
import gzip
import sys
import regex as re
import numpy as np
import math
import argparse
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools


def geteventstofilter(rMATStable):
	#Sometimes you want to remove events (say that are U2AF38 sensitive) from a set of events to then further analyze
	sensevents = []

	#Read U2AF38 table
	df = pd.read_table(rMATStable)

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

		if row['FDR'] < 0.05 and abs(row['IncLevelDifference']) >= 0.05:
			sensevents.append(eventID)

	print('{0} events were sensitive and will be filtered.'.format(len(sensevents)))

	return sensevents

def getSequences(rMATStable, genomefasta, sensevents):
	#Make (or blank) files
	for eventclass in ['moreincluded', 'lessincluded', 'control']:
		for seqregion in ['seq1', 'seq2', 'seq3', 'seq4']:
			outfh = open('{0}.{1}.fasta'.format(eventclass, seqregion), 'w')
			outfh.close()

	#Index fasta
	print('Indexing genome...')
	if os.path.basename(genomefasta).endswith('.gz'):
		seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta, 'rt'), 'fasta'))
	else:
		seq_dict = SeqIO.to_dict(SeqIO.parse(open(genomefasta, 'r'), 'fasta'))
	print('Done indexing!')

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
		sample1psis = row['IncLevel1'].split(',')
		sample2psis = row['IncLevel2'].split(',')
		if 'NA' in sample1psis or 'NA' in sample2psis:
			continue
		sample1psis = [float(psi) for psi in sample1psis]
		sample2psis = [float(psi) for psi in sample2psis]
		psis = sample1psis + sample2psis

		if row['FDR'] < 0.05 and row['IncLevelDifference'] >= 0.05 and (sensevents == None or (eventID not in sensevents)):
			eventclass = 'moreincluded'
		elif row['FDR'] < 0.05 and row['IncLevelDifference'] <= -0.05 and (sensevents == None or (eventID not in sensevents)):
			eventclass = 'lessincluded'
		elif row['FDR'] > 0.5 and max(psis) >= 0.15 and min(psis) <= 0.85 and (sensevents == None or (eventID not in sensevents)):
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

#Given a sequence and motif, where are the starts of motif instances within seq
def getmotifpos(seq, motif):
	starts = []
	seq = seq.upper()
	motif = motif.upper()
	pattern = r'{0}'.format(motif)
	for match in re.finditer(pattern, seq, overlapped = True):
		start = match.span()[0]
		starts.append(start)

	return starts


def searchfasta(fasta, motifs):
	seqsinfile = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqsinfile +=1

	seqcounter = 0
	fastamotifstartdict = {} #{seqname : [motif starts]}

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 500 == 0:
			print('Sequence {0} of {1}...'.format(seqcounter, seqsinfile))
		
		motifstarts_allmotifs = []
		for motif in motifs:
			motifstarts = getmotifpos(str(record.seq), motif)
			motifstarts_allmotifs += motifstarts
		#If there is no motif match, the motif start = 1000 (will be off the plot)
		if not motifstarts_allmotifs:
			motifstarts_allmotifs = [1000]
		
		
		#If this is seq2 or seq4, the positions need to be reversed.
		#i.e. they need to be counted from the right hand side of the sequence
		elif 'seq2' in fasta or 'seq4' in fasta:
			reversedmotifstarts_allmotifs = []
			for motifstart in motifstarts_allmotifs:
				reversedmotifstart = len(str(record.seq)) - motifstart
				reversedmotifstarts_allmotifs.append(reversedmotifstart)
			motifstarts = reversedmotifstarts_allmotifs
		

		fastamotifstartdict[str(record.id)] = motifstarts_allmotifs

	return fastamotifstartdict

def iteratefastas(motif):
	outfh = open('{0}results.txt'.format('gpattern'), 'w')
	outfh.close()

	fastas = []
	for eventtype in ['lessincluded', 'control', 'moreincluded']:
		for region in ['seq1', 'seq2', 'seq3', 'seq4']:
			fastas.append('{0}.{1}.fasta'.format(eventtype, region))

	#print fastas

	fastadfs = []
	for fasta in fastas:
		counter = 0
		dfs = []
		seqclass, region = fasta.split('.')[0], fasta.split('.')[1]
		fastamotifstartdict = searchfasta(fasta, motif) #{seqname : [motif starts]}
		#Add to DataFrame
		for seqname in fastamotifstartdict:
			for motifstart in fastamotifstartdict[seqname]:
				df = pd.DataFrame({'seqname' : [seqname], 'seqclass' : [seqclass], 'region' : [region], 'motifstart': [motifstart]}) #variables must be iterables for some reason
				dfs.append(df)

		fastadf = pd.concat(dfs, ignore_index = True)
		fastadfs.append(fastadf)
		
	bigdf = pd.concat(fastadfs, ignore_index = True)
	bigdf.to_csv('bigdftest.txt', sep = '\t', index = False)




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--filterrMATStable', type = str, help = 'Optional. rMATS table of events to filter from MATS table. These events will not be considered in further analyses.')
	parser.add_argument('--rMATStable', type = str, help = 'rMATS output table.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--motif', type = str, help = 'kmer or regular expression. If regular expression, surround with quotes and do not use 'r' at the beginning.')
	parser.add_argument('--scramblecontrol', action = 'store_true', help = 'If supplied, instead of looking for motif, look for all scrambles of the motif except the motif itself. Not useable if motif is a regular expression.')
	args = parser.parse_args()

	if args.filterrMATStable:
		sensevents = geteventstofilter(args.filterrMATStable)
	else:
		sensevents = None

	if args.scramblecontrol:
		motif = args.motif
	#Make all possible substrings of the motif of length k
	#Store these in subs
		k = 2
		l = len(motif)
		subs = []
		for i in range(l):
			for j in range(i, l):
				if (j+1) - i == k:
					subs.append(motif[i: j + 1])
		
		#all combinations, GC and CpG matched
		gccount = motif.count('G') + motif.count('C')
		cpgcount = motif.count('CG')
		allcombs = []
		for comb in itertools.product('ATCG', repeat = len(motif)):
			comb = ''.join(comb)
			combgccount = comb.count('G') + comb.count('C')
			combcpgcount = comb.count('CG')
			hassubmotif = False
			if combgccount == gccount and combcpgcount == cpgcount and comb != motif:
				for sub in subs:
					if sub in comb:
						hassubmotif = True
				if hassubmotif == False:
					allcombs.append(comb)

		motifs = allcombs



		
	elif not args.scramblecontrol:
		motifs = [args.motif]

	print('Analyzing content of these motifs: {0}'.format(motifs))

	getSequences(args.rMATStable, args.genomefasta, sensevents)
	iteratefastas(motifs)

#iteratefastas(r'(([AT]GG(.{0,7})[AT]GG(.{0,7})[AT]GG(.{0,7})[AT]GG))')


