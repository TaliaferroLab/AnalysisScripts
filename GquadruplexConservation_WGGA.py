#Using phylop scores, look at the conservation of G nucleotides contained with WGGA-type
#gquadruplex-forming motifs. Compare this to the conservation of other G nucleotides that are
#not in WGGA-type gquadruplex-forming motifs.  This is analogous to GquadruplexConservation.py, 
#except that this uses WGGA to look for G-quadruplexes while GquadruplexConservation.py uses
#RNAfold.

from Bio import SeqIO
import argparse
import regex
import os
import gffutils
import pysam
import numpy as np

def getallGpositions(seq):
	gpositions = []
	seq = seq.upper()
	for ind, nt in enumerate(seq):
		if nt == 'G':
			gpositions.append(ind)

	return gpositions

def getWGGAgpositions(seq):
	wggagpositions = [] #0-based coords of Gs in WGGA gquadruplex motifs
	seq = seq.upper().replace('U', 'T')
	motif = r'(([AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA(.{0,5})[AU]GGA))'

	#Get all overlapping instances of the motif with their start and end coords
	motifboundaries = [[m.start(), m.end()] for m in regex.finditer(motif, seq, overlapped = True)]
	for motifboundary in motifboundaries:
		motifstart = motifboundary[0]
		motifend = motifboundary[1]
		motifseq = seq[motifstart : motifend]
		for i in range(len(motifseq)):
			#Within each motif instance, get WGGA positions
			if motifseq[i:i+4] == 'AGGA' or motifseq[i:i+4] == 'TGGA':
				#The G positions have to be the middle to positions in WGGA
				gpositions = [motifstart + i + 1, motifstart + i + 2]
				for gposition in gpositions:
					if gposition not in wggagpositions:
						wggagpositions.append(gposition)

	return wggagpositions

def iteratefasta(fasta):
	gquadpositions = {} #{geneid : [wggagquadpositions, 0-based]}
	nongquadpositions = {} #{geneid : [nonwggagquadpositions, 0-based]}
	for record in SeqIO.parse(fasta, 'fasta'):
		genename = str(record.id).split('_')[0]
		seq = str(record.seq).upper().replace('U', 'T')
		nongquadgs = []
		allgpositions = getallGpositions(seq)
		wggagpositions = getWGGAgpositions(seq)
		for position in allgpositions:
			if position not in wggagpositions:
				nongquadgs.append(position)
		gquadpositions[genename] = wggagpositions
		nongquadpositions[genename] = nongquadgs

	return gquadpositions, nongquadpositions

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

	os.remove(db_fn)

	return gcoords

def getphylopscores(phylopbed, gcoords):
	scores = [] #all scores
	scoresd = {} #{geneid : [scores]}
	tbx = pysam.Tabixfile(phylopbed)
	nt = 0
	ntwithscores = 0

	for gene in gcoords:
		pcscores = [] #list of all pc scores for coords in this gene
		chrm = gcoords[gene][0]
		coords = gcoords[gene][2]
		for coord in coords:
			nt +=1
			for row in tbx.fetch(chrm, coord, coord + 1, parser = pysam.asTuple()):
				score = float(row[4])
				if score:
					ntwithscores +=1
					scores.append(score)
					pcscores.append(score)
		if pcscores:
			scoresd[gene] = pcscores

	print 'Interrogated {0} nucleotides. Found phylop scores for {1} ({2}%) of them.'.format(nt, ntwithscores, round((ntwithscores / float(nt)), 4) * 100)

	return scores, scoresd


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to consider.')
	parser.add_argument('--gff', type = str, help = 'Gff of regions.  Need it to get genomic coordinates.')
	parser.add_argument('--featurename', type = str, help = '3rd field of gff file of features you are interrogating. Usually UTR5, CDS, or UTR3.')
	parser.add_argument('--phylopbed', type = str, help = 'Phylop scores in gzipped, tabix indexed format.')
	args = parser.parse_args()

	gquadpositions, nongquadpositions = iteratefasta(args.fasta)
	gquadcoords = getcoords(gquadpositions, args.featurename, args.gff)
	nongquadcoords = getcoords(nongquadpositions, args.featurename, args.gff)
	scores, scoresd = getphylopscores(args.phylopbed, gquadcoords)
	print np.mean(scores)
	scores, scoresd = getphylopscores(args.phylopbed, nongquadcoords)
	print np.mean(scores)



