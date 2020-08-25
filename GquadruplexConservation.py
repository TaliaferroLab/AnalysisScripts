#For some reason on coyote, the fastest sqlite for use with gffutils is 3.7.13

import pysam
from Bio import SeqIO
import sys
import gffutils
import os
import pysam
import argparse

def getquadgs(gquadoutfile):
	#First get the coords of G's that are predicted to be in Gquads from gquadout.txt file (from UTRfold_Gquadruplex.py)
	gquadpositions = {} #{geneid : [positions in sequence]}
	with open(gquadoutfile, 'r') as gqf:
		for line in gqf:
			line = line.strip().split('\t')
			if line[0] == 'seqname':
				continue
			#gene ids look like ENSMUSG00000020097.14_ENSMUSG00000020097.14
			geneid = line[0].split('_')[0]
			positions = line[4]
			if positions != 'none':
				positions = [int(pos) for pos in positions.split(',')] #these positions are 0-based relative to the start of the sequence
			elif positions == 'none':
				positions = []
			gquadpositions[geneid] = positions

	return gquadpositions

def getnongquadgs(gquadoutfile, fasta):
	#Given a gquadoutfile and a fasta of the corresponding sequences, get the positions of any G that was not predicted to be in a quadruplex
	
	#Get positions of all G's in fasta
	allgpos = {} #{geneid: [all g positions, 0-based]}
	for record in SeqIO.parse(fasta, 'fasta'):
		geneid = str(record.id).split('_')[0]
		seq = str(record.seq)
		gpos = [pos for pos, nt in enumerate(seq) if nt == 'G']
		allgpos[geneid] = gpos

	nongquadpositions = {}
	with open(gquadoutfile, 'r') as gqf:
		for line in gqf:
			line = line.strip().split('\t')
			if line[0] == 'seqname':
				continue
			geneid = line[0].split('_')[0]
			positions = line[4]
			if positions != 'none':
				gquadpositions = [int(pos) for pos in positions.split(',')] #these positions are 0-based relative to the start of the sequence
			elif positions == 'none':
				gquadpositions = []
			#Remove anything from allgpos that is present in positions
			allg = allgpos[geneid]
			for gquadposition in gquadpositions:
				allg.remove(gquadposition)
			nongquadpositions[geneid] = allg

	return nongquadpositions


def getcoords(gpositions, feature, gff):
	#Given a list of gpositions within a fasta, get their genome coordinates
	#The gff is one that corresponds to the fasta, usually a longestXXX.gff3.
	#Feature is the the featuretype (3rd field in gff) that you want to go after, usually UTR5, CDS, or UTR3

	#Get all genome coords for a feature that are exonic
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	features = db.features_of_type(feature)
	featurecounter = 0

	exoniccoords = {} # {geneid : [chrm, strand, [list of all positions that are exonic, e.g. 400, 401, 402, 403, 410, 411, 412, etc.]]}
	#These will end up being 1-based gff-style coords
	for feature in features:
		featurecounter +=1
		if featurecounter % 5000 == 0:
			print 'Feature {0}...'.format(featurecounter)
		geneid = feature.attributes['Parent'][0].replace('gene:', '')
		chrm = str(feature.chrom)
		strand = str(feature.strand)
		exonicpos = []
		for exon in db.children(feature, featuretype = 'exon', level = 1):
			exonicnts = range(exon.start, exon.end + 1)
			exonicpos = exonicpos + exonicnts
		if strand == '+':
			exonicpos = sorted(exonicpos)
		elif strand == '-':
			exonicpos = sorted(exonicpos)[::-1]
		exoniccoords[geneid] = [chrm, strand, exonicpos]

	#Now intersect g positions within the feature with exonic coords
	gcoords = {} #{geneid : [chrm, [coords]]}
	for gene in gpositions:
		gpos = gpositions[gene]
		chrm = exoniccoords[gene][0]
		strand = exoniccoords[gene][1]
		ecoords = exoniccoords[gene][2]
		gc = []
		for gp in gpos:
			if strand == '+':
				gc.append(ecoords[gp])
			#If it's on the minus strand we have to count backwards from the end
			elif strand == '-':
				gc.append(ecoords[gp])
		gcoords[gene] = [chrm, strand, sorted(gc)]

	return gcoords

def getphastscores(phastconsbed, gcoords):
	scores = [] #all scores
	scoresd = {} #{geneid : [scores]}
	tbx = pysam.Tabixfile(phastconsbed)
	nt = 0
	ntwithscores = 0

	for gene in gcoords:
		pcscores = [] #list of all pc scores for coords in this gene
		chrm = gcoords[gene][0]
		coords = gcoords[gene][2]
		for coord in coords:
			nt +=1
			for row in tbx.fetch(chrm, coord, coord + 1, parser = pysam.asTuple()):
				score = row[4]
				if score:
					ntwithscores +=1
					scores.append(score)
					pcscores.append(score)
		if pcscores:
			scoresd[gene] = pcscores

	print 'Interrogated {0} nucleotides. Found phastcons scores for {1} ({2}%) of them.'.format(nt, ntwithscores, round((ntwithscores / float(nt)), 4) * 100)

	return scores, scoresd


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gtype', choices =['quad', 'nonquad'], type = str, help = 'Interrogate Gs participating in quads or not participating?')
	parser.add_argument('--gquadout', type = str, help = 'File produced by UTRfold_gquadruplex.py')
	parser.add_argument('--fasta', type = str, help = 'Fasta file of regions. Only necessary if gtype == nonquad. Need it to find out where the other Gs are.')
	parser.add_argument('--gff', type = str, help = 'Gff of regions. Need it to get coordinates.')
	parser.add_argument('--featurename', type = str, help = '3rd field of gff file of features you are interrogating. Usually UTR5, CDS, or UTR3.')
	parser.add_argument('--sequenceclass', type = str, help = 'Sequence bin. Usually AllGenes, Unchanged, SD1, SD2, or SD3.')
	parser.add_argument('--phastconsbed', type = str, help = 'Phastcons scores with tabix index in same directory.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	if not os.path.isfile(os.path.abspath(args.outfile)):
		with open(args.outfile, 'w') as f:
			f.write(('\t').join(['gene', 'score', 'gquad', 'Region', 'Class']) + '\n')

	if args.gtype == 'quad':
		gpositions = getquadgs(args.gquadout)
		gcoords = getcoords(gpositions, args.featurename, args.gff)
		scores, scoresd = getphastscores(args.phastconsbed, gcoords)

		with open(args.outfile, 'a') as f:
			for gene in scoresd:
				for score in scoresd[gene]:
					f.write(('\t').join([gene, str(score), 'true', args.featurename, args.sequenceclass]) + '\n')

	elif args.gtype == 'nonquad':
		gpositions = getnongquadgs(args.gquadout, args.fasta)
		gcoords = getcoords(gpositions, args.featurename, args.gff)
		scores, scoresd = getphastscores(args.phastconsbed, gcoords)

		with open(args.outfile, 'a') as f:
			for gene in scoresd:
				for score in scoresd[gene]:
					f.write(('\t').join([gene, str(score), 'false', args.featurename, args.sequenceclass]) + '\n')




