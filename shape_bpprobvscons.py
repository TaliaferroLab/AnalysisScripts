#Given an output of bp probabilities (from spats and RNAstructure) and a set of phastcons scores for the genome,
#output a matched list of bp probs and phastcons scores

import tabix
import pysam
import sys

def getbpprobsandcoords(reactivities, oligogff):
	#Probably use cutadapt2528_combined_pseudoT_fixed_normalized_reactivities_bpprobs.out
	bpprobs = {} # {oligoname : [{five_prime_offset : bpprob}, [chrm, start, stop, strand]]}
	fh = open(reactivities, 'r')
	for line in fh:
		line = line.strip().split('\t')
		if line[0] == 'sequence':
			continue
		oligoname = line[0]
		five_prime_offset = int(line[2])
		bpprob = float(line[10])
		if five_prime_offset >= 10 and five_prime_offset <= 90:
			if bpprobs.has_key(oligoname) == False:
				bpprobs[oligoname] = []
				bpprobs[oligoname].append({})
			bpprobs[oligoname][0][five_prime_offset] = bpprob

	fh.close()
	

	fh = open(oligogff, 'r')
	for line in fh:
		line = line.strip().split('\t')
		chrm = line[0]
		start = int(line[3])
		stop = int(line[4])
		strand = line[6]
		name = line[8][3:]
		name = name.split('|')
		name_noENS = ('|').join([name[0], name[1], name[2], name[4], name[5], name[6], name[7], name[8]])
		#if name_noENS in bpprobs:
			#bpprobs[name_noENS].append([chrm, start, stop, strand])
		name = ('|').join(name)
		if name in bpprobs and 'intron' in name and name.split('|')[0] == 'mammal':
			bpprobs[name].append([chrm, start, stop, strand])

	fh.close()
	return bpprobs

def getphastcons(phastconsbed, bpprobdict):
	phastconsdict = {} # {oligoname : {five_prime_offset : [bpprob, phastconsscore]}}
	#bpprobdict = {} # {oligoname : [{five_prime_offset : bpprob}, [chrm, start, stop, strand]]}
	phastconstabix = pysam.Tabixfile(phastconsbed)
	for oligo in bpprobdict:
		if len(bpprobdict[oligo]) == 1: #if there was no match for it in the gff (so it has not start/stop coords), skip it
			continue
		chrm, start, stop, strand = bpprobdict[oligo][1][0], bpprobdict[oligo][1][1], bpprobdict[oligo][1][2], bpprobdict[oligo][1][3]
		if strand == '+':
			try:
				for bed in phastconstabix.fetch(chrm, start-1, stop-1, parser = pysam.asBed()):
					if not (bed.start >= start and bed.end <= stop):
						continue
					phastconsscore = float(bed.name)
					fpo = int((bed.start - start) + 4) #the first nt of the gff region is five_prime_offset 4
					#print chrm, start, stop, bed.start, bed.end, fpo
					if phastconsdict.has_key(oligo) == False:
						phastconsdict[oligo] = {}
					if phastconsdict[oligo].has_key(fpo) == False and fpo >=10 and fpo <= 90: #discount anything near a splice site
						phastconsdict[oligo][fpo] = []
						bpprob = bpprobdict[oligo][0][fpo]
						phastconsdict[oligo][fpo] = [bpprob, phastconsscore]
			except ValueError:
				print 'WARNING: problem with {0}.'.format(oligo)

		elif strand == '-':
			try:
				for bed in phastconstabix.fetch(chrm, start-1, stop-1, parser = pysam.asBed()):
					if not (bed.start >= start and bed.end <= stop):
						continue
					phastconsscore = float(bed.name)
					fpo = int((stop - bed.start) + 3)
					#print chrm, start, stop, bed.start, bed.end, fpo
					if phastconsdict.has_key(oligo) == False:
						phastconsdict[oligo] = {}
					if phastconsdict[oligo].has_key(fpo) == False and fpo >= 10 and fpo <= 90: #discount anything near a splice site
						phastconsdict[oligo][fpo] = []
						bpprob = bpprobdict[oligo][0][fpo]
						phastconsdict[oligo][fpo] = [bpprob, phastconsscore]
			except ValueError:
				print 'WARNING: problem with {0}.'.format(oligo)

	return phastconsdict




bpprobdict = getbpprobsandcoords(sys.argv[1], sys.argv[2])
phastconsdict = getphastcons(sys.argv[3], bpprobdict)

print 'bpprob' + '\t' + 'phastcons'
for oligo in phastconsdict:
	for fpo in phastconsdict[oligo]:
		bpprob = phastconsdict[oligo][fpo][0]
		phastconsscore = phastconsdict[oligo][fpo][1]
		print str(bpprob) + '\t' + str(phastconsscore)

