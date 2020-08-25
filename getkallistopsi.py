import gffutils
import os
import sys
import re
from Bio import SeqIO
import itertools
import pickle
import gzip
import numpy as np
import random

def getpositionfactors(gff):
	genecount = 0
	txends = {} #{ENSMUSG : [strand, [list of distinct transcript end coords]]}
	posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Get number of distinct transcript ends for each gene
	genes = db.features_of_type('gene')
	for gene in genes:
		genename = str(gene.id).replace('gene:', '')
		ends = []
		if gene.strand == '+':
			for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'end'):
				if transcript.end not in ends:
					ends.append(transcript.end)
		elif gene.strand == '-':
			for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'start', reverse = True):
				if transcript.start not in ends:
					ends.append(transcript.start)

		if ends: #Sometimes there are no 'transcripts' for a gene, like with pseudogenes, etc.
			txends[genename] = [gene.strand, ends]

	#Sort transcript end coords
	s_txends = {} #{ENSMUSG : [sorted (most upstream to most downstream) tx end coords]}
	for gene in txends:
		strand = txends[gene][0]
		coords = txends[gene][1]
		if strand == '+':
			sortedcoords = sorted(coords)
		elif strand == '-':
			sortedcoords = sorted(coords, reverse = True)
		s_txends[gene] = sortedcoords


	#Figure out postion scaling factor for each transcript (position / (number of total positions - 1)) (m / (n - 1))
	genes = db.features_of_type('gene')
	for gene in genes:
		genecount +=1
		genename = str(gene.id).replace('gene:', '')
		if genename not in s_txends or len(s_txends[genename]) == 1:
			continue
		n = len(s_txends[genename])
		possibleends = s_txends[genename]
		posfactors[genename] = {}
		for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'end'):
			txname = str(transcript.id).replace('transcript:', '').split('.')[0]
			if gene.strand == '+':
				m = possibleends.index(transcript.end)
			elif gene.strand == '-':
				m = possibleends.index(transcript.start)
			posfactor = m / float(n - 1)
			posfactors[genename][txname] = posfactor

	print posfactors
	return posfactors

#Need to get UTR regions that are distinct to each isoform
#Start with gff of 3' UTR regions
def getdistinctregions(gff, genomefasta):
	distinctregions = {} #{geneid : {transcriptid(s) : [3UTR number, distinctUTRseq]}}
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		distinctseqs = {} #{transcriptid(s) : [pAsite counter (may be different than number of UTRs because not all UTRs are represented here, distinctUTRseq]}
		seenseqs = []
		utrcounter = 0
		mostdownstreamcoord = 0 #The most downstream coordinate of any UTR we've seen so far for this gene.
		geneid = str(gene.id).replace('gene:', '').split('.')[0]
		if gene.strand == '+':
			for UTR3 in db.children(gene, featuretype = 'UTR3', level = 1, order_by = 'end'):
				distinctseq = ''
				UTRid = str(UTR3.id).replace('UTR3:', '').split('.')[0]

				#If this is the first UTR for this gene
				if utrcounter == 0:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'start'):
						exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].upper()
						distinctseq += exonseq
					mostdownstreamcoord = UTR3.end
					utrcounter +=1
					distinctseqs[UTRid] = [utrcounter, str(distinctseq)]
				elif utrcounter >= 1:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'start'):
						#If this exon is somehow contained within the last one (should not be possible), skip it
						if exon.end <= mostdownstreamcoord:
							pass
						elif exon.end > mostdownstreamcoord:
							if exon.start < mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[mostdownstreamcoord:exon.end].upper()
								distinctseq += exonseq
							elif exon.start >= mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start - 1:exon.end].upper()
								distinctseq += exonseq
					
					mostdownstreamcoord = UTR3.end

					#Only going to call something a new polyA site if it's at least 50 nt away from the previous one
					#As a proxy for this, it must have at least 50 nt of "distinct" sequence
					if len(str(distinctseq)) >= 50:
						utrcounter +=1
						distinctseqs[UTRid] = [utrcounter, str(distinctseq)]

		elif gene.strand == '-':
			for UTR3 in db.children(gene, featuretype = 'UTR3', level = 1, order_by = 'start', reverse = True):
				distinctseq = ''
				UTRid = str(UTR3.id).replace('UTR3:', '').split('.')[0]

				#If this is the first UTR for this gene
				if utrcounter == 0:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'end', reverse = True):
						exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].reverse_complement().upper()
						#Must prepend instead of append this time
						distinctseq = distinctseq + exonseq
					mostdownstreamcoord = UTR3.start
					utrcounter +=1
					distinctseqs[UTRid] = [utrcounter, str(distinctseq)]
				elif utrcounter >= 1:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'end', reverse = True):
						#If this exon is somehow contained within the last one (should not be possible), skip it
						if exon.start >= mostdownstreamcoord:
							continue
						elif exon.start < mostdownstreamcoord:
							if exon.end > mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start-1:mostdownstreamcoord-1].reverse_complement().upper()
								distinctseq = distinctseq + exonseq
							elif exon.start <= mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].reverse_complement().upper()
								distinctseq = distinctseq + exonseq


					mostdownstreamcoord = UTR3.start
					if len(str(distinctseq)) >= 50:
						utrcounter +=1
						distinctseqs[UTRid] = [utrcounter, str(distinctseq)]

		distinctregions[geneid] = distinctseqs

	return distinctregions



#How many 3' UTRs does each gene have?
#Start with a gff of all 3' UTRs
def UTR3spergene(gff):
	utrs = {} #{geneid: [[UTRstart, UTRstop], [UTRstart, UTRstop]]}
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')
	for gene in genes:
		geneid = str(gene.id).replace('gene:', '')
		if geneid not in utrs:
			utrs[geneid] = []
		for utr in db.children(gene, level = 1, featuretype = 'UTR3'):
			if [utr.start, utr.end] not in utrs[geneid]:
				utrs[geneid].append([utr.start, utr.end])

	return utrs



#Given a sequence, return its g quadruplex density
def getgquaddens(seq, motif):
	seqlen = len(seq)
	#motifmatches = re.findall(motif, seq)
	#dens = len(motifmatches) / float(seqlen)
	dens = seq.count(motif) / float(seqlen)

	return dens

#Make all possible "versions" of the gquad motif
def makeallmotifs():
	motifs = []
	bases = ['A', 'T', 'G', 'C']
	#all3mers = [''.join(x) for x in itertools.product(bases, repeat = 4)]
	#for threemer in all3mers:
		#motif = r'(?=({0}(.{{0,7}}){0}(.{{0,7}}){0}(.{{0,7}}){0}))'.format(threemer)
		#motifs.append([motif, threemer])
	motifs = [''.join(x) for x in itertools.product(bases, repeat = 7)]

	return motifs


#Given a fasta (of all 3' UTRs?), get the gquad dens and then multiply that by its posfactor
def getscores(fasta, utrs, posfactors, motifs, outfile):
	with open(outfile, 'a') as f:
		f.write(('\t').join(['motif', 'gene', 'score']) + '\n')

	for m in motifs:
		motif = m[0]
		threemer = m[1]
		print motif, threemer
		scores = {} #{genename : [sum of gquaddens for all UTRs, [UTR1score, UTR2score]]}
		for record in SeqIO.parse(fasta, 'fasta'):
			txname = str(record.id).split('_')[0]
			genename = str(record.id).split('_')[1]
		
			#Only consider genes with >1 UTR
			if len(utrs[genename]) == 1:
				continue

			seq = str(record.seq)
			gquaddens = getgquaddens(seq, motif)
			posfactor = posfactors[txname][1]
			score = gquaddens * posfactor
			if genename not in scores:
				scores[genename] = [0, []]
			scores[genename][0] += gquaddens
			scores[genename][1].append(score)

		with open(outfile, 'a') as f:
			for gene in scores:
				totalscore = sum(scores[gene][1])
				totaldens = scores[gene][0]
				try:
					finalscore = totalscore / totaldens
				except ZeroDivisionError:
					finalscore = '0'
				f.write(('\t').join([threemer, gene, str(finalscore)]) + '\n')

def distinctregionscores(distinctregions, motifs, outfile):
	#distinctregions = {} #{geneid : {transcriptid(s) : [3UTR number, distinctUTRseq]}}
	with open(outfile, 'w') as f:
		f.write(('\t').join(['motif', 'gene', 'score']) + '\n')

	scores = {} #{motif : [scores]}
	for m in motifs:
		#motif = m[0]
		#threemer = m[1]
		motif = m
		threemer = m
		print motif, threemer
		for gene in distinctregions:
			totaldens = 0
			totalscore = 0
			#Only want genes with multiple distinct regions
			totalUTRs = len(distinctregions[gene])
			if totalUTRs == 1:
				continue
			for transcript in distinctregions[gene]:
				posfactor = (distinctregions[gene][transcript][0] - 1) / (float(totalUTRs) - 1)
				seq = distinctregions[gene][transcript][1]
				dens = getgquaddens(seq, motif)
				score = dens * posfactor
				totaldens += dens
				totalscore += score
			try:
				finalscore = totalscore / totaldens
				if motif not in scores:
					scores[motif] = [finalscore]
				else:
					scores[motif].append(finalscore)
			except ZeroDivisionError:
				finalscore = 'NA'
			with open(outfile, 'a') as f:
				f.write(('\t').join([threemer, gene, str(finalscore)]) + '\n')

	return scores

def getzscores(scores, zscoreout):
	meanscores = {} #{motif : meanscore across all genes}
	zscores = {} #{motif : zscore}
	for motif in scores:
		meanscore = np.mean(scores[motif])
		meanscores[motif] = meanscore
	motifcounter = 0
	for motif in meanscores:
		motifcounter +=1
		if motifcounter % 100 == 0:
			print 'Shuffling motif {0} of {1}...'.format(motifcounter, len(meanscores))
		shuffled = []
		shuffledmotifscores = []
		tries = 0
		cpgcount = motif.count('CG')
		l = list(motif)
		while len(shuffledmotifscores) <= 50 and tries <= 1000:
			tries +=1
			random.shuffle(l)
			shuffledmotif = ''.join(l)
			shuffledcpgcount = shuffledmotif.count('CG')
			if shuffledmotif != motif and shuffledmotif not in shuffled and cpgcount == shuffledcpgcount:
				shuffled.append(shuffledmotif)
				shuffledmotifscore = meanscores[shuffledmotif]
				shuffledmotifscores.append(shuffledmotifscore)

		if len(shuffled) == 0:
			zscore = 0
		else:
			shuffledsd = np.std(shuffledmotifscores)
			shuffledmean = np.mean(shuffledmotifscores)
			zscore = (meanscores[motif] - shuffledmean) / float(shuffledsd)
		zscores[motif] = zscore

	with open(zscoreout, 'a') as f:
		f.write(('\t').join(['motif', 'zscore']) + '\n')
		for motif in zscores:
			f.write(('\t').join([motif, str(zscores[motif])]) + '\n')











distinctregions = getdistinctregions(sys.argv[1], sys.argv[2])
#posfactors = getpositionfactors(sys.argv[1])
#utrs = UTR3spergene(sys.argv[2])
motifs = makeallmotifs()
#getscores(sys.argv[3], utrs, posfactors, motifs, sys.argv[4])
scores = distinctregionscores(distinctregions, motifs, sys.argv[3])
getzscores(scores, sys.argv[4])
