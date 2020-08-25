#Given a genome annotation (e.g. ensGene.mm9.gff3) and a list of genes you are interested in, retrieve 5' UTR coords (and/or sequences) for those genes.
#Defining 5' UTR regions is done as follows: The genome annotation has 5' UTR regions explicitly annotated (as five_prime_UTR).  There may be more than one 
#of these per gene, and they may be overlapping.  Define the 5' UTR region for a gene as beginning with the leftmost start coord of any five_prime_UTR entry
#for that gene and ending with the rightmost end coord of any five_prime_UTR entry. In short, get the most extreme boundaries for the five_prime_UTR entires 
#for that gene.

#Then, consider all exons for that gene, and take exonic regions that overlap or fall within the 5' UTR region defined above. Merge overlapping exons within
#this region to make longer exons. Finally

import gffutils
import sys
import os
from itertools import groupby
from operator import itemgetter
import argparse
from Bio import SeqIO

def ens2shortname(ens2shortfile):
	ens2short = {} #{ensembl gene id : shortname}
	infh = open(ens2shortfile, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			ens2short[line[0]] = line[1]

	infh.close()
	return ens2short

def getgenesofinterest(listofgenes):
	genesofinterest = []
	infh = open(listofgenes, 'r')
	for line in infh:
		line = line.strip()
		genesofinterest.append(line)
	infh.close()

	return genesofinterest

#Given a list of exon coordinates, merge overlapping exons into one long exon
#Derived from a SO post.
def mergeexons(exons):
	#exons = [[start, stop], [start, stop], [start, stop]]
	sortedexons = sorted(exons, key = lambda exon: exon[0])
	mergedexons = []

	for exon in sortedexons:
		if not mergedexons:
			mergedexons.append(exon)
		else:
			topexon = mergedexons.pop()
			if topexon[1] >= exon[0]:
				newexon = [topexon[0], max(topexon[1], exon[1])]
				mergedexons.append(newexon)
			else:
				mergedexons.append(topexon)
				mergedexons.append(exon)

	return mergedexons

#Given two intervals, return the overlapping integers between them
#Both exons are lists of coords [start, stop]
def getoverlap(exonA, exonB):
	A = range(exonA[0], exonA[1] + 1)
	B = range(exonB[0], exonB[1] + 1)
	A = set(A)
	overlap = list(A.intersection(B))

	return overlap


def longest5UTR(gff, ens2short, genesofinterest, outfile):
	analyzedgenes = 0 #counter
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'

	print 'Indexing gff...'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	print 'Done indexing!'

	db = gffutils.FeatureDB(db_fn)
	genes = db.features_of_type('gene')

	for gene in genes:
		UTR5exons = [] #[[start, end], [start, end]] for all exons that overlap extreme UTR5 boundaries
		CDSexons = [] # CDS regions
		filteredUTR5exons = [] #UTR5exon regions with any part that overlaps a CDS region removed
		lowestUTRstart = None
		highestUTRend = None
		
		#Only care about the genes that are of interest
		if gene.id not in ens2short:
			continue
		else:
			shortname = ens2short[gene.id]
			if shortname not in genesofinterest:
				continue
		
		#If there is no annotated 5' UTR region for this gene or if it has 0 coding exons, skip it
		five_prime_UTR_count = len(list(db.children(gene, featuretype = 'five_prime_UTR')))
		CDSexoncount = len(list(db.children(gene, featuretype = 'CDS')))
		if five_prime_UTR_count == 0 or CDSexoncount == 0:
			continue

		for UTR5 in db.children(gene, featuretype = 'five_prime_UTR'):
			#Get the most extreme coordinates for 5' UTR boundaries
			if not lowestUTRstart or UTR5.start < lowestUTRstart:
				lowestUTRstart = UTR5.start
			if not highestUTRend or UTR5.end > highestUTRend:
				highestUTRend = UTR5.end


		#Get all exon coordinates that lie in between the extreme 5' UTR boundaries

		for exon in db.children(gene, featuretype = 'exon'):
			#If the exon overhangs the extreme UTR start but not the extreme UTR end, the coords are from the UTR start to the exon end
			if exon.start < lowestUTRstart and exon.end <= highestUTRend and exon.end > lowestUTRstart:
				UTR5exons.append([lowestUTRstart, exon.end])
			#If the exon overhangs the extreme UTR end but not the extreme UTR start, the coords are from the exon start to the UTR end
			elif exon.start >= lowestUTRstart and exon.start < highestUTRend and exon.end > highestUTRend:
				UTR5exons.append([exon.start, highestUTRend])
			#If the exon boundaries are contained within the extreme UTR boundaries, the coords are from the exon start to the exon end
			elif exon.start >= lowestUTRstart and exon.end <= highestUTRend:
				UTR5exons.append([exon.start, exon.end])


		#Get all CDS exon boundaries. Doesn't matter if they overlap or are repeated as we are going to remove any UTR5 exon portion if it overlaps any CDS exon portion.
		for CDS in db.children(gene, featuretype = 'CDS'):
			CDSexons.append([CDS.start, CDS.end])

		#Merge overlapping exons
		UTR5exons = mergeexons(UTR5exons)

		for UTR5exon in UTR5exons:
			overlap = [] #list of coords that are overlaps between UTR5exons and CDSexons
			for CDSexon in CDSexons:
				#Get the overlapping positions between this UTR exon and this CDS exon
				overlap += getoverlap(UTR5exon, CDSexon)
				overlap = sorted(overlap)
			#If there is no overlap between CDS exons and 5UTRexons
			if not overlap:
				if UTR5exon not in filteredUTR5exons:
					filteredUTR5exons.append(UTR5exon)
				continue
			else:
				exon = set(range(UTR5exon[0], UTR5exon[1]))
				overlap = set(overlap)
				filtered = sorted(list(exon - overlap))
				#It could be true that an overlap is entirely contained within the UTR5exon. If this is the case, simply taking the first and last values of "filtered"
				#will just return the original UTR5exon coords with the undesired CDS overlap still there.  To get around this, look for groups of consecutive integers in filtered.
				#Then just take the first and last values of each group of consecutive integers.  Most of the time, there will just be one group of consecutive integers because the
				#CDS will not be entirely contained within the UTR5exon.
				if filtered:
					for k, g in groupby(enumerate(filtered), lambda (i, x): i - x):
						f = map(itemgetter(1), g)
						if len(f) > 1:
							UTR5exon = [f[0], f[-1]]
							if UTR5exon not in filteredUTR5exons:
								filteredUTR5exons.append(UTR5exon)


		if len(filteredUTR5exons) >= 1:
			analyzedgenes +=1
			with open(outfile, 'a') as f:
				#These filteredUTR5exons may not actually be exons. They are exonic sequence, but have been merged together and have had CDS overlaps removed.
				#We are going to write each one out individually because I don't want to join them together, creating artficial kmers at the junctions.
				filteredUTR5exons = sorted(filteredUTR5exons, key = itemgetter(0))
				#Write "gene" level line, which in this case is the extreme coordiates of the 5' UTR
				f.write(('\t').join([gene.chrom, 'UTR5s', 'UTR5', str(filteredUTR5exons[0][0]), str(filteredUTR5exons[-1][1]), '.', gene.strand, '.', 'Name=' + gene.id + ';ID=' + shortname]) + '\n')
				for ind, exon in enumerate(filteredUTR5exons):
					f.write(('\t').join([gene.chrom, 'UTR5s', 'exon', str(exon[0]), str(exon[1]), '.', gene.strand, '.', 'Name=' + gene.id + '.exon' + str(ind) + ';Parent=' + shortname + ';ID=' + shortname + '.exon' + str(ind)]) + '\n')

	print 'Found 5\' UTRs for {0} of {1} genes.'.format(analyzedgenes, len(genesofinterest))

	#os.remove(db_fn)

#Get the sequences of the UTR regions. Do NOT join "exons" together. Instead, since exons are "merged" exon coordinates that have had overlapping coding regions removed,
#report each "exon" as its own fasta entry.

def gff2seq(gff, genomeseq, output):
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomeseq, 'fasta'))
	gfffh = open(gff, 'r')
	for line in gfffh:
		line = line.strip().split('\t')
		if line[2] == 'exon':
			chrm = line[0]
			start = int(line[3])
			stop = int(line[4])
			strand = line[6]
			name = line[8].split(';ID=')[-1]
			if 'random' in chrm:
				continue

			if strand == '+':
				seq = seq_dict[chrm].seq[start - 1 : stop].upper()
			elif strand == '-':
				seq = seq_dict[chrm].seq[start - 1 : stop].reverse_complement().upper()

			with open(output, 'a') as f:
				f.write('>' + name + '\n' + str(seq) + '\n')

	gfffh.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--genesofinterest', type = str, help = 'File containing list of genes you are interested in. Names should be gene short names.', required = True)
	parser.add_argument('--ens2short', type = str, help = 'File containing ensembl gene names and their corresponding gene short names. Usually Ensembl2GeneName_Mouse.txt.', required = True)
	parser.add_argument('--gff', type = str, help = 'Gff file of genome annotation with five_prime_UTR and CDS exons explicitly annotated. Usually ensGene.mm9.gff3.', required = True)
	parser.add_argument('--outputgff', type = str, help = 'Output gff of 5\' UTR exons.', required = True)
	parser.add_argument('--givemeseqs', action = 'store_true', help = 'Do you also want a fasta file of the UTRs?')
	parser.add_argument('--genomeseq', type = str, help = 'Genome sequence in fasta format. Required if using givemeseqs.')
	parser.add_argument('--outputfasta', type = str, help = 'Output file for UTR seqs. Required if using givemeseqs.')
	args = parser.parse_args()

	ens2short = ens2shortname(args.ens2short)
	genesofinterest = getgenesofinterest(args.genesofinterest)
	
	#Reset any file with the output name to a blank file.
	fh = open(args.outputgff, 'w')
	fh.close()
	longest5UTR(args.gff, ens2short, genesofinterest, args.outputgff)
	
	if args.givemeseqs:
		if not args.genomeseq or not args.outputfasta:
			print 'You must supply both the genome sequence in fasta format and an output file for your sequences!'
			sys.exit()
		fh = open(args.outputfasta, 'w')
		fh.close()
		gff2seq(args.outputgff, args.genomeseq, args.outputfasta)


