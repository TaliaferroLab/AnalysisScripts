#Like GquadruplexConservation.py but instead of calculating conservation using PhastCons or PhyloP scores,
#look to see if syntenic nt in human is G.
#If Library not loaded: libssl.1.0.0.dylib: export DYLD_LIBRARY_PATH=/usr/local/opt/openssl/lib/ ; source ~/.bash_profile
#This works best with MySQL 5.6, rather than the newest version 5.7.
#5.7 doesn't seem to work with mysql-python 1.2.5
#MySQL installed with homebrew and mysql-python installed with anaconda

from cogent.db.ensembl import Genome, Compara
import os
import gffutils
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import gzip

def indexgenome(genomefasta):
	print 'Indexing genome...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print '{0} chromosomes indexed.'.format(len(seq_dict))

	return seq_dict


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

	os.remove(db_fn)

	return gcoords


#A large percentage of the time, the returned alignment is of the wrong strand (i.e. the reverse complement).
#We need to see if this is true by comparing to genome sequence.
#If it is true, then take the reverse complements of both the mouse and human alignments
def checkstrand(mousealign, humanalign, chrm, coords, seq_dict):
	possiblemouseseq = mousealign.replace('-', '')
	genomeseq = seq_dict[chrm].seq[min(coords)-1:max(coords)]
	if possiblemouseseq == genomeseq:
		#then everything is cool
		return mousealign, humanalign
	elif possiblemouseseq != genomeseq:
		mouseseq = Seq(mousealign, generic_dna)
		mouseseq = str(mouseseq.reverse_complement()).replace('-', '')
		if mouseseq == genomeseq and humanalign != 'None':
			#conveniently, the rev comp of '-' is '-'
			mousealign = str(Seq(mousealign, generic_dna).reverse_complement())
			humanalign = str(Seq(humanalign, generic_dna).reverse_complement())
			return mousealign, humanalign
		elif mouseseq == genomeseq and humanalign == 'None':
			mousealign = str(Seq(mousealign, generic_dna).reverse_complement())
			humanalign = 'None'
			return mousealign, humanalign
		elif mouseseq != genomeseq:
			print 'ERROR: sequences aren\'t matching even after reverse complementing!'
			return None, None



#Need to have a way of handling gaps to see if sequences are the same
#Can figure out which non '-' character we are interested in based on coords
#Then get index of that character in string with '-'s
#Then see if that index is the same character in human sequence
def checkifmatch(mouseseq, humanseq, coords):
	indices = [coord - min(coords) for coord in coords]
	gappedindices = [] #gquad coords of gapped seq string, e.g. GACT---ATGCA---ACACG
	matches = [] #e.g. [True, True, False, False, True] in the order they are in coords (which should be sorted)
	nongapcounter = -1
	
	#First check to see if there is any aligned human seq. If there isn't, it is 'None'.
	if humanseq == 'None':
		return [False] * len(coords)

	for ind, nt in enumerate(mouseseq):
		if nt != '-':
			nongapcounter +=1
			if nongapcounter in indices:
				gappedindices.append(ind)

	for ind in gappedindices:
		if mouseseq[ind] == humanseq[ind]:
			match = True
		elif mouseseq[ind] != humanseq[ind]:
			match = False
		matches.append(match)

	return matches


def getsynteny(gcoords):
	print 'Looking for synteny...'
	genecounter = 0
	compara = Compara(['mouse', 'human'], Release = 83, account = None)
	gcoords_matches = {} #just like gcoords but with a list of matches {gene : [chrm, strand, [coords], [matches]]}
	for gene in gcoords:
		genecounter +=1
		if genecounter % 100 == 0:
			print 'Analyzing gene {0} of {1}...'.format(genecounter, len(gcoords))
		info = gcoords[gene]
		chrm = info[0].replace('chr', '') #ensembl style chromosome names
		#Strand is not used by compara.getSyntenicRegions
		strand = info[1]
		coords = info[2]
		#If there were no g coordinates for this gene
		if not coords:
			continue
		#Get alignment for entire length of region using min(coords) and max(coords) to minimize number of calls.
		#Parse out the results later.
		for synt_region in compara.getSyntenicRegions(Species = 'mouse', CoordName = chrm, Start = min(coords), End = max(coords), Strand = strand, ensembl_coord = True, align_method = 'PECAN', align_clade = '23 amniota vertebrates Pecan'):
			membs = synt_region.Members
			if len(membs) != 2:
				#Not really sure how this can happen.  Usually if there is no human alignment, it returns None
				continue
			mouse = membs[0]
			human = membs[1]
			mouseseq = str(mouse.AlignedSeq)
			humanseq = str(human.AlignedSeq)
			mouseseq, humanseq = checkstrand(mouseseq, humanseq, 'chr' + chrm, coords, seq_dict)
			if mouseseq and humanseq:
				matches = checkifmatch(mouseseq, humanseq, coords)
				info.append(matches)
				#Not every gene makes it this far.  Some have no syntentic regions I think.
				gcoords_matches[gene] = info
			else:
				print 'Sequence problems with {0}.  Alignment and genome sequence did not match.'.format(gene)
	
	return gcoords_matches

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gtype', choices =['quad', 'nonquad'], type = str, help = 'Interrogate Gs participating in quads or not participating?')
	parser.add_argument('--gquadout', type = str, help = 'File produced by UTRfold_gquadruplex.py')
	parser.add_argument('--fasta', type = str, help = 'Fasta file of regions. Only necessary if gtype == nonquad. Need it to find out where the other Gs are.')
	parser.add_argument('--gff', type = str, help = 'Gff of regions. Need it to get coordinates.')
	parser.add_argument('--featurename', type = str, help = '3rd field of gff file of features you are interrogating. Usually UTR5, CDS, or UTR3.')
	parser.add_argument('--sequenceclass', type = str, help = 'Sequence bin. Usually AllGenes, Unchanged, SD1, SD2, or SD3.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	if not os.path.isfile(os.path.abspath(args.outfile)):
		with open(args.outfile, 'w') as f:
			f.write(('\t').join(['gene', 'chrom', 'strand', 'class', 'region', 'coord', 'gquad', 'in_human']) + '\n')

	seq_dict = indexgenome(args.genomefasta)


	if args.gtype == 'quad':
		gpositions = getquadgs(args.gquadout)
		gcoords = getcoords(gpositions, args.featurename, args.gff)
		gcoords_matches = getsynteny(gcoords)

	elif args.gtype == 'nonquad':
		gpositions = getnongquadgs(args.gquadout, args.fasta)
		gcoords = getcoords(gpositions, args.featurename, args.gff)
		gcoords_matches = getsynteny(gcoords)

	with open(args.outfile, 'a') as f:
		for gene in gcoords_matches:
			info = gcoords_matches[gene]
			chrm = info[0]
			strand = info[1]
			coords = info[2]
			matches = info[3]
			for ind, coord in enumerate(coords):
				match = matches[ind]
				f.write(('\t').join([gene, chrm, strand, args.sequenceclass, args.featurename, str(coord), args.gtype, str(match)]) + '\n')



