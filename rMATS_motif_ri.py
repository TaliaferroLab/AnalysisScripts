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

def geteventstofilter(rMATStable):
	#Sometimes you want to remove events (say that are U2AF38 sensitive) from a set of events to then further analyze
	sensevents = []

	#Read U2AF38 table
	df = pd.read_table(rMATStable)

	for index, row in df.iterrows():
		chrm = row['chr']
		strand = row['strand']
		riexonstart = row['riExonStart_0base']
		riexonend = row['riExonEnd']
		upstreamexonstart = row['upstreamES']
		upstreamexonend = row['upstreamEE']
		downstreamexonstart = row['downstreamES']
		downstreamexonend = row['downstreamEE']
		eventID = chrm + ':' + strand + ':' + str(riexonstart) + '-' + str(riexonend) + ':' + str(upstreamexonstart) + '-' + str(upstreamexonend) + ':' + str(downstreamexonstart) + '-' + str(downstreamexonend)
		FDR = row['FDR']
		dpsi = row['IncLevelDifference']

		if row['FDR'] < 0.05 and abs(row['IncLevelDifference']) >= 0.05:
			sensevents.append(eventID)

	print '{0} events were sensitive and will be filtered.'.format(len(sensevents))

	return sensevents

def getSequences(rMATStable, genomefasta, sensevents):
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

	#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| (this is on the plus strand)
	#|||||||||||||||||||--------------------------------------------------------------------------------------------|||||||||||||||||||
	#__________________  _________________                                                         _________________ _________________
	#  seq1 -150 to 0  seq2 0 to 150                                                               seq3 -150 to 0    seq4 0 to 150
	#            <stretch1>     <stretch2>                                                             <stretch3>          <stretch4>

	for index, row in df.iterrows():
		chrm = row['chr']
		strand = row['strand']
		riexonstart = row['riExonStart_0base']
		riexonend = row['riExonEnd']
		upstreamexonstart = row['upstreamES']
		upstreamexonend = row['upstreamEE']
		downstreamexonstart = row['downstreamES']
		downstreamexonend = row['downstreamEE']
		eventID = chrm + ':' + strand + ':' + str(riexonstart) + '-' + str(riexonend) + ':' + str(upstreamexonstart) + '-' + str(upstreamexonend) + ':' + str(downstreamexonstart) + '-' + str(downstreamexonend)
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
			stretch1 = min(150, upstreamexonend - upstreamexonstart)
			stretch2 = min(150, downstreamexonstart - upstreamexonend)
			stretch3 = min(150, downstreamexonstart - upstreamexonend)
			stretch4 = min(150, downstreamexonend - downstreamexonstart)

		elif strand == '-':
			stretch1 = min(150, downstreamexonend - downstreamexonstart)
			stretch2 = min(150, downstreamexonstart - upstreamexonend)
			stretch3 = min(150, downstreamexonstart - upstreamexonend)
			stretch4 = min(150, upstreamexonend - upstreamexonstart)


		if strand == '+':
			seq1 = seq_dict[chrm].seq[upstreamexonend - stretch1 : upstreamexonend]
			seq2 = seq_dict[chrm].seq[upstreamexonend : upstreamexonend + stretch2]
			seq3 = seq_dict[chrm].seq[downstreamexonstart - stretch3 : downstreamexonstart]
			seq4 = seq_dict[chrm].seq[downstreamexonstart : downstreamexonstart + stretch4]

		elif strand == '-':
			seq1 = seq_dict[chrm].seq[downstreamexonstart : downstreamexonstart + stretch1].reverse_complement()
			seq2 = seq_dict[chrm].seq[downstreamexonstart - stretch2 : downstreamexonstart].reverse_complement()
			seq3 = seq_dict[chrm].seq[upstreamexonend : upstreamexonend + stretch3].reverse_complement()
			seq4 = seq_dict[chrm].seq[upstreamexonend - stretch4 : upstreamexonend].reverse_complement()

		for fn, seqregion in zip(['seq1', 'seq2', 'seq3', 'seq4'], [str(seq1).upper(), str(seq2).upper(), str(seq3).upper(), str(seq4).upper()]):
			with open('{0}.{1}.fasta'.format(eventclass, fn), 'a') as outfh:
				outfh.write('>' + eventID + '\n' + seqregion + '\n')

#Given a sequence and motif, find how many nucleotides in the seq are covered by the motif
def getmotifpos(seq, motif):
	coverednt = []
	seq = seq.upper()
	motif = motif.upper()
	pattern = r'{0}'.format(motif)
	for match in re.finditer(pattern, seq, overlapped = True):
		coveredbymatch = range(match.span()[0], match.span()[1])
		coverednt += coveredbymatch
	coverednt = list(set(coverednt))
	coverednt = sorted(coverednt)

	return coverednt

def searchfasta(fasta, windowsize, slidesize, motif):
	seqsinfile = 0
	for record in SeqIO.parse(fasta, 'fasta'):
		seqsinfile +=1

	seqcounter = 0
	fastamotifmatchdict = {} #{seqname : {0-based window center : fraction of nt covered by motif in this window}}
	fastamotifcovereddict = {} #{seqname : {0-based position : <0 if not covered by motif in any window, 1 if covered by motif in at least one window>}}

	for record in SeqIO.parse(fasta, 'fasta'):
		seqcounter +=1
		if seqcounter % 500 == 0:
			print 'Sequence {0} of {1}...'.format(seqcounter, seqsinfile)
		motifdict = {} #{0-basedposition window center : fraction of nt covered by motif in this window}

		#Create and initialize motifcovereddict
		motifcovereddict = {} #{0-based position : <0 if not covered by motif in any window, 1 if covered by motif in at least one window>}
		for i in range(200):
			motifcovereddict[record.id] = 0
		seq = str(record.seq)


		windowstart = 0
		while windowstart + windowsize <= len(seq):
			windowseq = seq[windowstart : windowstart + windowsize]
			windowcenter = windowstart + (len(windowseq) / 2)
			coverednt = getmotifpos(windowseq, motif)
			motifdict[windowcenter] = len(coverednt) / float(len(windowseq))

			#Adjust position of motif locations by window position
			coverednt = [pos + windowstart for pos in coverednt]
			#if this position is covered by a motif, it becomes 1 in motifcovereddict
			for pos in coverednt:
				motifcovereddict[pos] = 1

			windowstart += slidesize

		fastamotifmatchdict[record.id] = motifdict
		fastamotifcovereddict[record.id] = motifcovereddict

	#Write gquadpos per seq out to file
	outfile = os.path.splitext(fasta)[0] + '.{0}pos.txt'.format('gpattern')
	with open(outfile, 'w') as outfh:
		for seqid in fastamotifcovereddict:
			motifpos = []
			for pos in sorted(fastamotifcovereddict[seqid]):
				if fastamotifcovereddict[seqid][pos] == 1:
					motifpos.append(pos)
			motifpos = sorted(motifpos)
			motifpos = [str(pos) for pos in motifpos]
			motifpos = (',').join(motifpos)
			if not motifpos:
				motifpos = 'None'
			outfh.write(seqid + '\t' + motifpos + '\n')

	#print fastamotifmatchdict
	return fastamotifmatchdict

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
		fastamotifmatchdict = searchfasta(fasta, 40, 1, motif) #{seqname : {0-based window center : fraction of nt covered by motif in this window}}
		#Add to DataFrame
		for seqname in fastamotifmatchdict:
			counter +=1
			if counter % 50 == 0:
				print counter
			for windowcenter in sorted(fastamotifmatchdict[seqname]):
				df = pd.DataFrame({'seqname' : [seqname], 'seqclass' : [seqclass], 'region' : [region], 'pos': [windowcenter], 'values' : [fastamotifmatchdict[seqname][windowcenter]]}) #variables must be iterables for some reason
				dfs.append(df)

		fastadf = pd.concat(dfs, ignore_index = True)
		fastadfs.append(fastadf)
		fastadf.to_csv('{0}.{1}.df.txt'.format(seqclass, region), sep = '\t', index = False)
		
	bigdf = pd.concat(fastadfs, ignore_index = True)
	bigdf.to_csv('bigdf.txt', sep = '\t', index = False)

def plotbigdf():
	from ggplot import *
	bigdf = pd.read_table('bigdf.txt')
	#Group by position, region, and seqclass
	bigdf = bigdf.groupby(['pos', 'region', 'seqclass'], as_index = False)
	#Get the mean at each position/region/seqclass
	means = bigdf.agg(['mean', 'count', 'std'])
	#Reset index (which was made during grouping)
	means = means.reset_index()
	means.columns = ['pos', 'region', 'seqclass', 'mean', 'count', 'std']

	#Now plot
	p = ggplot(means, aes(x = 'pos', y = 'mean', color = 'seqclass')) + geom_line() + facet_wrap('region')

	#Save
	p.save(filename = './means.pdf')

	#Statistical tests
	bigdf = pd.read_table('bigdf.txt')
	#Get unique values
	posvalues = bigdf.pos.unique()
	regions = bigdf.region.unique()
	seqclasses = bigdf.seqclass.unique()
	statdfs = [] 

	for pos in posvalues: #1, 2, 3, etc.
		for region in regions: #seq1, seq2, etc.
			valuesdict = {}
			for seqclass in seqclasses: #lessincluded, moreincluded, etc.
				values = bigdf[(bigdf['pos'] == pos) & (bigdf['region'] == region) & (bigdf['seqclass'] == seqclass)]
				valuesdict[seqclass] = values
			try:
				lessincludedp = mannwhitneyu(valuesdict['lessincluded']['values'], valuesdict['control']['values'])[1]
			except ValueError: #if all values are identical
				lessincludedp = 1.0
			try:
				moreincludedp = mannwhitneyu(valuesdict['moreincluded']['values'], valuesdict['control']['values'])[1]
			except ValueError:
				moreincludedp = 1.0


			df = pd.DataFrame({'pos' : [pos], 'region' : [region], 'sample' : ['moreincludedp'], 'pval' : [moreincludedp]})
			statdfs.append(df)
			df = pd.DataFrame({'pos' : [pos], 'region' : [region], 'sample' : ['lessincludedp'], 'pval' : [lessincludedp]})
			statdfs.append(df)

	#Combine statdfs
	statdf = pd.concat(statdfs, ignore_index = True)

	#Multiple hypothesis correct
	pvals = statdf['pval']
	correctedpvals = multipletests(pvals, method = 'fdr_bh')[1]
	logcorrectedpvals = [np.log10(pval) * -1 for pval in correctedpvals]
	df = pd.DataFrame({'correctedpval' : correctedpvals})
	df2 = pd.DataFrame({'logcorrectedpval' : logcorrectedpvals})
	statdf = pd.concat([statdf, df, df2], axis = 1) #cbind


	#plot
	p = ggplot(statdf, aes(x = 'pos', y = 'np.log10(pval) * -1', color = 'sample')) + geom_line() + facet_wrap('region')
	p.save(filename = './pvalues.pdf')





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--filterrMATStable', type = str, help = 'Optional. rMATS table of events to filter from MATS table. These events will not be considered in further analyses.')
	parser.add_argument('--rMATStable', type = str, help = 'rMATS output table.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--motif', type = str, help = 'kmer or regular expression. If regular expression, surround with quotes and do not use 'r' at the beginning.')
	args = parser.parse_args()

	if args.filterrMATStable:
		sensevents = geteventstofilter(args.filterrMATStable)
	else:
		sensevents = None

	getSequences(args.rMATStable, args.genomefasta, sensevents)
	iteratefastas(args.motif)
	plotbigdf()