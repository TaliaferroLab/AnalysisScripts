import pandas as pd
import os
import sys
from Bio import SeqIO
import gzip


def preptable(table):
	#Want to do a little bit of prep on a rMATS output table.  Add in some columns, split others.
	#Read in table
	df = pd.read_table(table)

	
	#Split replicate psi values so they have their own columns
	#This splits the column delimited replicate values so that each one has its own column
	#df = df.join(df['IncLevel1'].str.split(',', expand = True).rename(columns = {0: 'S1.Rep1', 1: 'S1.Rep2', 2: 'S1.Rep3'}))
	#df = df.join(df['IncLevel2'].str.split(',', expand = True).rename(columns = {0: 'S2.Rep1', 1: 'S2.Rep2', 2: 'S2.Rep3'}))
	

	#IDs are not stable from run to run so we have to make a new one
	eventIDs = []
	for index, row in df.iterrows():
		chrm = row['chr']
		strand = row['strand']
		exonstart = str(row['exonStart_0base'])
		exonend = str(row['exonEnd'])
		upstreamexonstart = str(row['upstreamES'])
		upstreamexonend = str(row['upstreamEE'])
		downstreamexonstart = str(row['downstreamES'])
		downstreamexonend = str(row['downstreamEE'])
		eventID = chrm + ':' + strand + ':' + upstreamexonstart + '-' + upstreamexonend + ':' + exonstart + '-' + exonend + ':' + downstreamexonstart + '-' + downstreamexonend
		eventIDs.append(eventID)
	df['eventID'] = eventIDs

	
	#Add in the total number of reads (inc isoform + exc isoform) for each replicate
	#Then add a column that says whether or not each replicate passes a read filter threshold
	s1readcounts = []
	s2readcounts = []
	readcountpass = []
	for index, row in df.iterrows():
		s1Inc = row['IJC_SAMPLE_1'].split(',')
		s1Sk = row['SJC_SAMPLE_1'].split(',')
		s2Inc = row['IJC_SAMPLE_2'].split(',')
		s2Sk = row['SJC_SAMPLE_2'].split(',')
	
		#Turn them into integers
		s1Inc = [int(x) for x in s1Inc]
		s1Sk = [int(x) for x in s1Sk]
		s2Inc = [int(x) for x in s2Inc]
		s2Sk = [int(x) for x in s2Sk]

		s1reads = [] #[rep1readcount, rep2readcount, rep3readcount] readcount = inclusion reads + exclusion reads
		s2reads = [] #[rep1readcount, rep2readcount, rep3readcount]
		for i in range(len(s1Inc)):
			s1reads.append(s1Inc[i] + s1Sk[i])
		for i in range(len(s2Inc)):
			s2reads.append(s2Inc[i] + s2Sk[i])

		#We need at least 20 reads in every replicate of both samples for this quantification to be believable
		if all(x >= 20 for x in s1reads) and all(x >= 20 for x in s2reads):
			readcountpass.append('pass')
		else:
			readcountpass.append('fail')

		#Turn them into strings
		s1reads = [str(x) for x in s1reads]
		s2reads = [str(x) for x in s2reads]

		s1readcounts.append((',').join(s1reads))
		s2readcounts.append((',').join(s2reads))

	df['s1reads'] = s1readcounts
	df['s2reads'] = s2readcounts
	df['readcountpass'] = readcountpass
	

	outfile = os.path.splitext(table)[0] + '.prepped.txt'
	df.to_csv(outfile, sep = '\t', index = False)

#Given a rMATS output table, get relevant sequences surrounding skipped exons
def getSequences(table, genomefasta):
	#Index genome
	print 'Indexing genome...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	seqs = {}
	for eventclass in ['moreincluded', 'lessincluded', 'control']:
		seqs[eventclass] = {} #{exonid : {upstream : seq, downstream : seq, exon: seq}}

	#Go through table
	df = pd.read_table(table)

	for index, row in df.iterrows():
		#if row['FDR'] < 0.05 and row['readcountpass'] == 'pass' and row['IncLevelDifference'] >= 0.05:
		if row['FDR'] < 0.05 and row['IncLevelDifference'] > 0:
			eventclass = 'moreincluded'
		#elif row['FDR'] < 0.05 and row['readcountpass'] == 'pass' and row['IncLevelDifference'] <= -0.05:
		elif row['FDR'] < 0.05 and row['IncLevelDifference'] < 0:
			eventclass = 'lessincluded'
		#elif row['FDR'] > 0.05 and row['readcountpass'] == 'pass' and abs(row['IncLevelDifference']) <= 0.05:
		elif row['FDR'] > 0.05:
			eventclass = 'control'
		else:
			continue

		ID = row['ID']
		geneid = row['GeneID']
		exonid = geneid + '_' + str(ID)
		chrm = row['chr']
		strand = row['strand']
		exonstart = row['exonStart_0base']
		exonend = row['exonEnd']

		if strand == '+':
			exonseq = seq_dict[chrm].seq[exonstart : exonend]
			upstreamseq = seq_dict[chrm].seq[exonstart - 100 : exonstart]
			downstreamseq = seq_dict[chrm].seq[exonend : exonend + 100]
		elif strand == '-':
			exonseq = seq_dict[chrm].seq[exonstart : exonend].reverse_complement()
			downstreamseq = seq_dict[chrm].seq[exonstart - 100 : exonstart].reverse_complement()
			upstreamseq = seq_dict[chrm].seq[exonend: exonend + 100].reverse_complement()

		#Add seq to dictionary
		seqs[eventclass][exonid] = {}
		seqs[eventclass][exonid]['upstream'] = str(upstreamseq).upper()
		seqs[eventclass][exonid]['downstream'] = str(downstreamseq).upper()
		seqs[eventclass][exonid]['exon'] = str(exonseq).upper()

	print 'Found sequences for {0} more included exons, {1} less included exons, and {2} control exons.'.format(len(seqs['moreincluded']), len(seqs['lessincluded']), len(seqs['control']))

	with open('Control_upstream.fasta', 'w') as upstreamfh, open('Control_exon.fasta', 'w') as exonfh, open('Control_downstream.fasta', 'w') as downstreamfh:
		for exonid in seqs['control']:
			upstreamfh.write('>' + exonid + '\n' + seqs['control'][exonid]['upstream'] + '\n')
			exonfh.write('>' + exonid + '\n' + seqs['control'][exonid]['exon'] + '\n')
			downstreamfh.write('>' + exonid + '\n' + seqs['control'][exonid]['downstream'] + '\n')

	with open('MoreIncluded_upstream.fasta', 'w') as upstreamfh, open('MoreIncluded_exon.fasta', 'w') as exonfh, open('MoreIncluded_downstream.fasta', 'w') as downstreamfh:
		for exonid in seqs['moreincluded']:
			upstreamfh.write('>' + exonid + '\n' + seqs['moreincluded'][exonid]['upstream'] + '\n')
			exonfh.write('>' + exonid + '\n' + seqs['moreincluded'][exonid]['exon'] + '\n')
			downstreamfh.write('>' + exonid + '\n' + seqs['moreincluded'][exonid]['downstream'] + '\n')

	with open('LessIncluded_upstream.fasta', 'w') as upstreamfh, open('LessIncluded_exon.fasta', 'w') as exonfh, open('LessIncluded_downstream.fasta', 'w') as downstreamfh:
		for exonid in seqs['lessincluded']:
			upstreamfh.write('>' + exonid + '\n' + seqs['lessincluded'][exonid]['upstream'] + '\n')
			exonfh.write('>' + exonid + '\n' + seqs['lessincluded'][exonid]['exon'] + '\n')
			downstreamfh.write('>' + exonid + '\n' + seqs['lessincluded'][exonid]['downstream'] + '\n')






preptable(sys.argv[1])
#getSequences(sys.argv[1], sys.argv[2])