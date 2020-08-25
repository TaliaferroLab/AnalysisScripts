from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import gzip
import numpy as np
import argparse
import gffutils
import os
import re
import random
import subprocess
import subsamplegff


def sam2gff(sam, gffoutfile):
	#sam is the mapping of the original sequences from junjie
	#gffoutfile is a gff file we will make of those alignments
	#Each sequence was 61 nt long
	#Some cross splice junctions

	analyzedseqs = 0
	seqsthatpass = 0

	with open(sam, 'r') as infile, open(gffoutfile, 'w') as gffout:
		for line in infile:
			analyzedseqs +=1
			line = line.strip().split('\t')
			ID = line[0]
			if line[1] == '0':
				strand = '+'
			elif line[1] == '16':
				strand = '-'
			else:
				continue
			chrm = line[2]
			start = int(line[3])
			
			cigar = line[5]
			if 'N' in cigar:
				spliced = True
			elif 'N' not in cigar:
				spliced = False
			#We aren't going to deal with anything that has an insertion or a deletion or any clipping
			if 'I' in cigar or 'D' in cigar or 'S' in cigar or 'H' in cigar:
				continue

			seqsthatpass +=1
			seq = Seq(line[9], generic_dna)
			if strand == '+':
				seq = str(seq)
			elif strand == '-':
				seq = str(seq.reverse_complement())

			#If it's not spliced, writing the gff is easy.
			if spliced == False:
				regionattributes = 'ID={0};Seq={1}'.format(ID, seq)
				exonattributes = 'ID={0}.exon1;Parent={1};Seq={2}'.format(ID, ID, seq)
				gffout.write(('\t').join([chrm, 'junjieG4', 'region', str(start), str(start + 61), '.', strand, '.', regionattributes]) + '\n')
				gffout.write(('\t').join([chrm, 'junjieG4', 'exon', str(start), str(start + 61), '.', strand, '.', exonattributes]) + '\n')

			#If it is spliced, it's slightly more complicated
			elif spliced == True:
				exon1length = int(cigar.split('M')[0])
				intronlength = int(cigar.split('M')[1].split('N')[0])
				exon2length = int(cigar.split('N')[1].split('M')[0])
				exon1start = start
				exon1end = start + exon1length - 1
				exon2start = start + exon1length + intronlength
				exon2end = start + exon1length + intronlength + exon2length - 1
				exon1seq = seq[:exon1length]
				exon2seq = seq[exon2length*-1:]
				
				regionattributes = 'ID={0};Seq={1}'.format(ID, seq)
				exon1attributes = 'ID={0}.exon1;Parent={1};Seq={2}'.format(ID, ID, exon1seq)
				exon2attributes = 'ID={0}.exon2;Parent={1};Seq={2}'.format(ID, ID, exon2seq)

				gffout.write(('\t').join([chrm, 'junjieG4', 'region', str(start), str(exon2end), '.', strand, '.', regionattributes]) + '\n')
				gffout.write(('\t').join([chrm, 'junjieG4', 'exon', str(exon1start), str(exon1end), '.', strand, '.', exon1attributes]) + '\n')
				gffout.write(('\t').join([chrm, 'junjieG4', 'exon', str(exon2start), str(exon2end), '.', strand, '.', exon2attributes]) + '\n')

	print 'Analyzed {0} seqs. {1} of these passed filters and were written to a gff.'.format(analyzedseqs, seqsthatpass)

def getGratio(seqs):
	#Get the ratio of Gs to Cs in a fasta

	ratios = []
	for seq in seqs:
		seq = seq.upper()
		g = seq.count('G')
		c = seq.count('C')
		if g == 0 or c == 0:
			continue
		else:
			ratio = float(g) / (float(g) + float(c))
			ratios.append(ratio)

	return ratios

def getcGcCscore(seq):
	#Given a sequence, calculate its cG/cC score.
	#First, get the longest consecutive substring of Gs and of Cs
	#http://stackoverflow.com/questions/18776238/count-the-number-of-max-consecutive-as-from-a-string-python-3

	seq = seq.upper()
	if seq.count('G') == 0:
		maxG = 0
	else:
		maxG = max(len(s) for s in re.findall(r'G+', seq))

	if seq.count('C') == 0:
		maxC = 0
	else:
		maxC = max(len(s) for s in re.findall(r'C+', seq))

	longestrun = max(maxG, maxC)

	cGscore = 0
	cCscore = 0
	#First get the cG score
	for i in range(1, longestrun + 1):
		searchstring = 'G' * i
		matches = re.findall(r'(?=({0}))'.format(searchstring), seq)
		score = len(matches) * i
		cGscore += score

	#Now the cC score
	for i in range(1, longestrun + 1):
		searchstring = 'C' * i
		matches = re.findall(r'(?=({0}))'.format(searchstring), seq)
		score = len(matches) * i
		cCscore += score

	if cCscore == 0:
		cGcCscore = cGscore
	else:
		cGcCscore = cGscore / float(cCscore)

	return cGcCscore

def getWGGAdens(seq):
	#Given a sequence, return its WGGA-(N0-6)-WGGA-(N0-6)-WGGA-(N0-6)-WGGA dens
	seqlen = len(seq)
	motifmatches = re.findall(r'(?=([AUT]GGA(.{0,5})[AUT]GGA(.{0,5})[AUT]GGA(.{0,5})[AUT]GGA))', seq)
	dens = len(motifmatches) / float(seqlen)
	return dens

def fold_a_fasta_gquad(fasta):
	#Given a fasta of gquad (or control) regions, fold each region and count the number of G's participating in gquads

	gquads = 0 #cumulative number of g's participating in gquads in this fasta
	energies = [] #list of folding energies for seqs in this fasta

	command = 'RNAfold -g < {0}'.format(fasta)
	job = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE,
		stderr = subprocess.PIPE, stdin = subprocess.PIPE)
	output = job.communicate()[0][1:] #take first value of tuple and remove first character ('>')
	perseq = output.split('>')
	for seq in perseq:
		structure = seq.split('\n')[2].split(' ')[0]
		energy = seq.split('\n')[2].split(' (')[1].replace(' ', '')[:-1]
		gquads += structure.count('+')
		energies.append(energy)

	for f in os.listdir('.'):
		if f.endswith('.ps'):
			os.remove(f)

	return gquads, energies

def makecontrolfasta(gquadfasta, alltxfasta):
	gquadtxs = []
	sampledseqs = []
	#First get all tx names that have a gquad in them (usually from mESCGquad.fasta)
	for record in SeqIO.parse(gquadfasta, 'fasta'):
		ID = str(record.id).split('_')
		ID = ('_').join([ID[0], ID[1]])
		gquadtxs.append(ID)

	#Now randomly sample 61 bp regions from these transcripts
	for record in SeqIO.parse(gzip.open(alltxfasta), 'fasta'):
		ID = str(record.id)
		if ID in gquadtxs:
			seq = str(record.seq)
			seqlen = len(seq)
			#We can't have a start position that is within 61nt of the end
			possiblestarts = range(0, seqlen - 61)
			start = random.choice(possiblestarts)
			sampledseq = seq[start : start + 61]
			sampledseqs.append(sampledseq.upper())

	return sampledseqs

def comparegratios(gquadfasta, alltxfasta, outfile):
	gratiodict = {} #{'sample' : [gratios], 'control' : [gratios]}
	sampleseqs = []

	for record in SeqIO.parse(gquadfasta, 'fasta'):
		seq = str(record.seq.upper())
		sampleseqs.append(seq)

	sampleratios = getGratio(sampleseqs)
	gratiodict['sample'] = sampleratios
	gratiodict['control'] = []

	#Make 100 different control seqs for each transcript in gquadfasta
	i = 0
	while i < 100:
		if (i + 1) % 10 == 0:
			print i + 1
		controlseqs = makecontrolfasta(gquadfasta, alltxfasta)
		gratios = getGratio(controlseqs)
		gratiodict['control'] += gratios
		i += 1

	with open(outfile, 'w') as f:
		f.write('sample' + '\t' + 'gratio' + '\n')
		for samp in gratiodict:
			for gratio in gratiodict[samp]:
				f.write(samp + '\t' + str(gratio) + '\n')

def comparecGcCscore(gquadfasta, alltxfasta, outfile):
	cGcCdict = {} #{'sample' : [scores], 'control' : [scores]}
	samplescores = []

	for record in SeqIO.parse(gquadfasta, 'fasta'):
		seq = str(record.seq.upper())
		score = getcGcCscore(seq)
		samplescores.append(score)

	cGcCdict['sample'] = samplescores
	cGcCdict['control'] = []

	#Make 100 different control seqs for each transcript in gquadfasta
	i = 0
	while i < 100:
		if (i + 1) % 10 == 0:
			print i + 1
		controlscores = []
		controlseqs = makecontrolfasta(gquadfasta, alltxfasta)
		for seq in controlseqs:
			score = getcGcCscore(seq)
			controlscores.append(score)
		#Add this batch of control scores to the dictionary
		cGcCdict['control'] += controlscores
		i += 1

	with open(outfile, 'w') as f:
		f.write('sample' + '\t' + 'cGcCscore' + '\n')
		for samp in cGcCdict:
			for score in cGcCdict[samp]:
				f.write(samp + '\t' + str(score) + '\n')

def compareWGGAdens(gquadfasta, alltxfasta, outfile):
	WGGAdict = {} #{'sample' : [densities], 'control' : [densities]}
	sampledensities = []

	for record in SeqIO.parse(gquadfasta, 'fasta'):
		seq = str(record.seq.upper())
		density = getWGGAdens(seq)
		sampledensities.append(density)

	WGGAdict['sample'] = sampledensities
	WGGAdict['control'] = []

	#Make 100 different control seqs for each transcript in gquadfasta
	i = 0
	while i < 100:
	 	if (i + 1) % 10 == 0:
	 		print i + 1
	 	controldensities = []
	 	controlseqs = makecontrolfasta(gquadfasta, alltxfasta)
	 	for seq in controlseqs:
	 		density = getWGGAdens(seq)
	 		controldensities.append(density)
	 	#Add this batch of control densities to the dictionary
	 	WGGAdict['control'] += controldensities
	 	i += 1

	#Can we bootstrap error for what fraction of sample and control populations contain a WGGA?
	fracs = []
	for i in range(1000):
		densities = []
		for i in range(int(round(len(sampledensities) * 0.2))):
			sampleddensity = random.choice(sampledensities)
			densities.append(sampleddensity)

		#True for every nonzerodensity, False for every zero density
		nonzerodensities = [True if density > 0 else False for density in densities]
		frac = sum(nonzerodensities) / float(len(nonzerodensities))
		fracs.append(frac)
	print np.mean(fracs), np.std(fracs), np.std(fracs) / float(len(fracs)**0.5)

	controldensities = WGGAdict['control']
	fracs = []
	for i in range(1000):
		densities = []
		for i in range(int(round(len(sampledensities) * 0.2))):
			sampleddensity = random.choice(controldensities)
			densities.append(sampleddensity)

		#True for every nonzerodensity, False for every zero density
		nonzerodensities = [True if density > 0 else False for density in densities]
		frac = sum(nonzerodensities) / float(len(nonzerodensities))
		fracs.append(frac)
	print np.mean(fracs), np.std(fracs), np.std(fracs) / float(len(fracs)**0.5)




	with open(outfile, 'w') as f:
	 	f.write('sample' + '\t' + 'WGGAdens' + '\n')
	 	for samp in WGGAdict:
	 		for dens in WGGAdict[samp]:
	 			f.write(samp + '\t' + str(dens) + '\n')

def comparernafold(gquadfasta, alltxfasta, densityoutfile, energyoutfile):
	#Given a gquad fasta, fold all the seqs with RNAfold and count how many G's are participating in quads.
	#Then make a 100 control fastas (one at a time) using seqs from the same transcripts and do the same thing.
	#These seqs should already be short (~60 nt), so there is no need to chop them.
	foldingenergies = {} #{'g4gquads' : [list of folding energies of seqs], 'control' : [list of folding energies of seqs]}

	#First, get number of seqs in gquadfasta
	numberofseqs = 0
	for record in SeqIO.parse(gquadfasta, 'fasta'):
		numberofseqs +=1

	#Get number of G's in gquadfasta that are participating in quads
	print 'Folding gquad fasta...'
	g4quads, energies = fold_a_fasta_gquad(gquadfasta)
	g4quaddensity = g4quads / float(numberofseqs)
	foldingenergies['g4gquads'] = energies
	foldingenergies['control'] = []


	#Make a fasta from the same txs as those that gquad regions are from
	#Fold it.
	#Do this 100 times.
	controlgquaddensities = []
	for i in range(100):
		if (i + 1) % 10 == 0:
			print 'Folding control fasta {0}...'.format(i + 1)
		seqs = makecontrolfasta(gquadfasta, alltxfasta)
		with open('temp.fasta', 'w') as f:
			for i in range(len(seqs)):
				f.write('>' + 'ControlSeq{0}'.format(i+1) + '\n' + seqs[i] + '\n')
		controlgquads, energies = fold_a_fasta_gquad('temp.fasta')
		foldingenergies['control'] += energies
		controlgquaddensity = controlgquads / float(numberofseqs)
		controlgquaddensities.append(controlgquaddensity)

	with open(densityoutfile, 'w') as f:
		f.write('sample' + '\t' + 'G4gsPerseq' + '\n')
		f.write('gquadfasta' + '\t' + str(g4quaddensity) + '\n')
		for controldensity in controlgquaddensities:
			f.write('control' + '\t' + str(controldensity) + '\n')

	with open(energyoutfile, 'w') as f:
		f.write('sample' + '\t' + 'foldingenergy' + '\n')
		for sample in foldingenergies:
			for energy in foldingenergies[sample]:
				f.write(sample + '\t' + str(energy) + '\n')




def overlapgffs(g4gff, txregiongff, modelgff):
	#Make gff databases
	print 'Indexing gff...'
	gff_fn = g4gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	g4db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing gff...'
	gff_fn = txregiongff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True, merge_strategy = 'merge')
	txregiondb = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Subsample txregiongff so that its length looks like modelgff
	#If you want to verify that the subsampling worked, look at temp.txt
	subsampledfeatures = subsamplegff.subsampleLength(modelgff, txregiongff, 40, 'temp.txt')
	#os.remove('temp.txt')

	g4overlaps = {} #{txregionid : [list of g4regionids that overlap]}
	genecounter = 0
	overlapcounter = 0
	totaltxregionlength = 0
	txregiongenes = txregiondb.features_of_type('gene')
	for gene in txregiongenes:
		if str(gene.id) not in subsampledfeatures:
			continue
		txg4overlap = False
		genecounter +=1
		if genecounter % 500 == 0:
			print 'Analyzing gene {0} in txregiongff...'.format(genecounter)
		g4overlaps[str(gene.id)] = []
		#Get all the exons for this txregion
		exoncoords = [] #[[exon1start, exon1end], [exon2start, exon2end]]
		for exon in txregiondb.children(gene, featuretype = 'exon', level = 2):
			exoncoords.append([exon.start, exon.end])
			exonlength = exon.end - exon.start
			totaltxregionlength += exonlength

		g4regions = g4db.features_of_type('region')
		for g4region in g4regions:
			foundoverlap = False
			if g4region.chrom == gene.chrom and g4region.strand == gene.strand:
				for g4exon in g4db.children(g4region, featuretype = 'exon', level = 1):
					g4exonnt = range(g4exon.start, g4exon.end + 1)
					for txexon in exoncoords:
						txexonnt = range(txexon[0], txexon[1] + 1)
						overlap = list(set(g4exonnt).intersection(txexonnt))
						if len(overlap) > 0:
							g4overlaps[str(gene.id)].append(str(g4region.id))
							foundoverlap = True
							txg4overlap = True
							#overlapcounter += len(overlap)

							#One overlap per g4 region
							overlapcounter +=1
							break

					#If you found an overlap with this exon of the g4 region, that's enough
					#Don't look in the next exon of this g4 region to avoud doublecounting
					if foundoverlap == True:
						break




	#print g4overlaps
	#print overlapcounter / float(genecounter)
	#print overlapcounter, genecounter, totaltxregionlength, overlapcounter / float(totaltxregionlength)
	return overlapcounter / float(genecounter)


def overlapgffs_bootstrap(g4gff, txregiongff):
	txregionlengths = {} #{txregion : length}
	#Make gff databases
	print 'Indexing gff...'
	gff_fn = g4gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	g4db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing gff...'
	gff_fn = txregiongff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True, merge_strategy = 'merge')
	txregiondb = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	g4overlaps = {} #{txregionid : [nt that overlap g4 region]}
	genecounter = 0
	overlapcounter = 0
	totaltxregionlength = 0
	txregiongenes = txregiondb.features_of_type('gene')
	for gene in txregiongenes:
		txregionlength = 0
		txg4overlappingnt = 0
		txg4overlap = False
		genecounter +=1
		if genecounter % 500 == 0:
			print 'Analyzing gene {0} in txregiongff...'.format(genecounter)
		g4overlaps[str(gene.id)] = 0
		#Get all the exons for this txregion
		exoncoords = [] #[[exon1start, exon1end], [exon2start, exon2end]]
		for exon in txregiondb.children(gene, featuretype = 'exon', level = 2):
			exoncoords.append([exon.start, exon.end])
			exonlength = exon.end - exon.start
			txregionlength += exonlength
			totaltxregionlength += exonlength

		txregionlengths[str(gene.id)] = txregionlength

		g4regions = g4db.features_of_type('region')
		for g4region in g4regions:
			foundoverlap = False
			if g4region.chrom == gene.chrom and g4region.strand == gene.strand:
				for g4exon in g4db.children(g4region, featuretype = 'exon', level = 1):
					g4exonnt = range(g4exon.start, g4exon.end)
					for txexon in exoncoords:
						txexonnt = range(txexon[0], txexon[1])
						overlap = list(set(g4exonnt).intersection(txexonnt))
						if len(overlap) > 0:
							#Add number of g4 overlap nucleotides divided by txregion exon length to dictionary
							txg4overlappingnt += len(overlap)
							foundoverlap = True
							txg4overlap = True
							overlapcounter += len(overlap)
							break

					#If you found an overlap with this exon of the g4 region, that's enough
					#Don't look in the next exon of this g4 region to avoid doublecounting
					#if foundoverlap == True:
						#break

		g4overlaps[str(gene.id)] = txg4overlappingnt

	#print g4overlaps
	densities = [] #list of all the densities of subsampled tx regions
	i = 0
	while i < 100:
		sampledvalues = []
		sampledlengths = []
		numbertosubsample = int(round(len(g4overlaps) * 0.3))
		possibletxs = g4overlaps.keys()
		for j in range(numbertosubsample):
			sampledtx = random.choice(possibletxs)
			value = g4overlaps[sampledtx]
			txlen = txregionlengths[sampledtx]
			sampledvalues.append(value)
			sampledlengths.append(txlen)
		density = np.sum(sampledvalues) / float(np.sum(sampledlengths))
		densities.append(density)
		i +=1
	print np.mean(densities), (np.std(densities) / float(i)**0.5)
	print max(g4overlaps.values())
	#print np.mean(g4overlaps.values())
	print overlapcounter, totaltxregionlength, overlapcounter / float(totaltxregionlength)





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', choices = ['sam2gff', 'gratio', 'cGcC', 'WGGA', 'rnafold', 'overlapgffs', 'overlapbootstrap'], required = True)
	parser.add_argument('--sam', type = str, help = 'Sam file of G4 regions.')
	parser.add_argument('--gffout', type = str, help = 'Gff output of G4 regions.')
	parser.add_argument('--gquadfasta', type = str, help = 'Input fasta file.')
	parser.add_argument('--alltxfasta', type = str, help = 'All transcripts (usually refseq so that name matches with gquadfasta names).')
	parser.add_argument('--g4gff', type = str, help = 'Experimental G4 regions in gff format.')
	parser.add_argument('--txregiongff', type = str, help = 'Regions of interest in gff format.')
	parser.add_argument('--modelgff', type = str, help = 'We want to subsample txregiongff so that its lengths match the distribution of those in modelgff.')
	parser.add_argument('--region', type = str)
	parser.add_argument('--txclass', type = str)
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	parser.add_argument('--energyoutfile', type = str, help = 'Output file for folding energies.  Needed if mode is rnafold.')
	args = parser.parse_args()

	if args.mode == 'sam2gff':
		sam2gff(args.sam, args.gffout)

	elif args.mode == 'gratio':
		comparegratios(args.gquadfasta, args.alltxfasta, args.outfile)

	elif args.mode == 'cGcC':
		comparecGcCscore(args.gquadfasta, args.alltxfasta, args.outfile)

	elif args.mode == 'WGGA':
		compareWGGAdens(args.gquadfasta, args.alltxfasta, args.outfile)

	elif args.mode == 'rnafold':
		comparernafold(args.gquadfasta, args.alltxfasta, args.outfile, args.energyoutfile)

	elif args.mode == 'overlapgffs':
		densities = []
		for i in range(100):
			if (i + 1) % 10 == 0:
				print 'Subsample {0} of 100...'.format(i + 1)
			density = overlapgffs(args.g4gff, args.txregiongff, args.modelgff)
			densities.append(density)

		with open('G4overlapdensities_2foldcutoff.txt', 'a') as f:
			for density in densities:
				f.write(args.txclass + '\t' + args.region + '\t' + str(density) + '\n')

		print np.mean(densities), np.std(densities), np.std(densities) / float(len(densities) ** 0.5)

	elif args.mode == 'overlapbootstrap':
		overlapgffs_bootstrap(args.g4gff, args.txregiongff)




