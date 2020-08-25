#Given a tpm table from kallisto data, calculate "psi" values by weighting expression of transcripts based on the 3' end coordinate.
#See Book 3 pg 139 for a graphical depiction of this "psi".

import gffutils
import os
import sys
from Bio import SeqIO
import itertools
import gzip
import argparse



def getpositionfactors(gff):
	genecount = 0
	txends = {} #{ENSMUSG : [strand, [list of distinct transcript end coords]]}
	posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Get number of distinct transcript ends for each gene
	genes = db.features_of_type('gene')
	for gene in genes:
		genename = str(gene.id).replace('gene:', '').split('.')[0]
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
		genename = str(gene.id).replace('gene:', '').split('.')[0]
		if genename not in s_txends or len(s_txends[genename]) == 1: #any pseudogene or other that has no transcripts or len(s_txends[genename]) == 1:
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

	return posfactors

def gettpmtable(tpmtable):
	#like ~/Documents/MIT/Localizaiton/Fmr1Deletion/kallisto/Gencodecomp/tpmtable.txt
	tpms ={} #{ensgid : {enstid : {WTSoma : tpm, WTNeurite : tpm, KOSoma : tpm, KONeurite : tpm}}}
	with open(tpmtable, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'target_id':
				continue
			enstid = line[0]
			ensgid = line[1]
			WTsomaAtpm = float(line[4])
			WTsomaBtpm = float(line[5])
			WTsomaCtpm = float(line[6])
			WTneuriteAtpm = float(line[7])
			WTneuriteBtpm = float(line[8])
			WTneuriteCtpm = float(line[9])
			KOsomaAtpm = float(line[10])
			KOsomaBtpm = float(line[11])
			KOsomaCtpm = float(line[12])
			KOneuriteAtpm = float(line[13])
			KOneuriteBtpm = float(line[14])
			KOneuriteCtpm = float(line[15])
			if ensgid not in tpms:
				tpms[ensgid] = {}
			if enstid not in tpms[ensgid]:
				tpms[ensgid][enstid] = {}
			tpms[ensgid][enstid]['WTSomaA'] = WTsomaAtpm
			tpms[ensgid][enstid]['WTSomaB'] = WTsomaBtpm
			tpms[ensgid][enstid]['WTSomaC'] = WTsomaCtpm
			tpms[ensgid][enstid]['WTNeuriteA'] = WTneuriteAtpm
			tpms[ensgid][enstid]['WTNeuriteB'] = WTneuriteBtpm
			tpms[ensgid][enstid]['WTNeuriteC'] = WTneuriteCtpm
			tpms[ensgid][enstid]['KOSomaA'] = KOsomaAtpm
			tpms[ensgid][enstid]['KOSomaB'] = KOsomaBtpm
			tpms[ensgid][enstid]['KOSomaC'] = KOsomaCtpm
			tpms[ensgid][enstid]['KONeuriteA'] = KOneuriteAtpm
			tpms[ensgid][enstid]['KONeuriteB'] = KOneuriteBtpm
			tpms[ensgid][enstid]['KONeuriteC'] = KOneuriteCtpm

	txcount = 0
	for gene in tpms:
		for tx in tpms[gene]:
			txcount +=1

	print 'There are {0} genes and {1} transcripts in this tpmtable.'.format(len(tpms), txcount)
	return tpms

def gettotaltpms(tpms):
	#For each gene, get the sum of all tpms across transcripts for each sample
	totaltpms = {} #{ensgid : {WTSoma : tpm, WTNeurite : tpm, KOSoma : tpm, KONeurite : tpm}}
	for gene in tpms:
		WTsomaAtpm = 0
		WTsomaBtpm = 0
		WTsomaCtpm = 0
		WTneuriteAtpm = 0
		WTneuriteBtpm = 0
		WTneuriteCtpm = 0
		KOsomaAtpm = 0
		KOsomaBtpm = 0
		KOsomaCtpm = 0
		KOneuriteAtpm = 0
		KOneuriteBtpm = 0
		KOneuriteCtpm = 0
		for tx in tpms[gene]:
			WTsomaAtpm += tpms[gene][tx]['WTSomaA']
			WTsomaBtpm += tpms[gene][tx]['WTSomaB']
			WTsomaCtpm += tpms[gene][tx]['WTSomaC']
			WTneuriteAtpm += tpms[gene][tx]['WTNeuriteA']
			WTneuriteBtpm += tpms[gene][tx]['WTNeuriteB']
			WTneuriteCtpm += tpms[gene][tx]['WTNeuriteC']
			KOsomaAtpm += tpms[gene][tx]['KOSomaA']
			KOsomaBtpm += tpms[gene][tx]['KOSomaB']
			KOsomaCtpm += tpms[gene][tx]['KOSomaC']
			KOneuriteAtpm += tpms[gene][tx]['KONeuriteA']
			KOneuriteBtpm += tpms[gene][tx]['KONeuriteB']
			KOneuriteCtpm += tpms[gene][tx]['KONeuriteC']
		totaltpms[gene] = {}
		totaltpms[gene]['WTSomaA'] = WTsomaAtpm
		totaltpms[gene]['WTSomaB'] = WTsomaBtpm
		totaltpms[gene]['WTSomaC'] = WTsomaCtpm
		totaltpms[gene]['WTNeuriteA'] = WTneuriteAtpm
		totaltpms[gene]['WTNeuriteB'] = WTneuriteBtpm
		totaltpms[gene]['WTNeuriteC'] = WTneuriteCtpm
		totaltpms[gene]['KOSomaA'] = KOsomaAtpm
		totaltpms[gene]['KOSomaB'] = KOsomaBtpm
		totaltpms[gene]['KOSomaC'] = KOsomaCtpm
		totaltpms[gene]['KONeuriteA'] = KOneuriteAtpm
		totaltpms[gene]['KONeuriteB'] = KOneuriteBtpm
		totaltpms[gene]['KONeuriteC'] = KOneuriteCtpm

	return totaltpms


def getpsis(posfactors, tpms, totaltpms):
	psis = {} #{ensgid : {WTSoma : psi, WTNeurite : psi, KOSoma : psi, KONeurite : psi}}

	for compartment in ['WTSomaA', 'WTSomaB', 'WTSomaC', 'WTNeuriteA', 'WTNeuriteB', 'WTNeuriteC', 'KOSomaA', 'KOSomaB', 'KOSomaC', 'KONeuriteA', 'KONeuriteB', 'KONeuriteC']:
		for gene in tpms:
			weightedtpms = []
			totaltpm = totaltpms[gene][compartment]
			#Only care about genes in posfactors (which is ones with >1 3' end)
			if gene not in posfactors:
				continue
			for tx in tpms[gene]:
				posfactor = posfactors[gene][tx]
				tpm = tpms[gene][tx][compartment]
				weightedtpm = posfactor * tpm
				weightedtpms.append(weightedtpm)
			#Only want genes that are expressed in at least one compartment
			if totaltpm < 5:
				continue
			psi = sum(weightedtpms) / float(totaltpm)
			if gene not in psis:
				psis[gene] = {}
			psis[gene][compartment] = psi

	return psis


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Genome annotation from which to get 3\' UTR order.')
	parser.add_argument('--tpmtable', type = str, help = 'TPM estimates in table format.  See ~/Documents/MIT/Localizaiton/Fmr1Deletion/kallisto/Gencodecomp/tpmtable.txt for example.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	posfactors = getpositionfactors(args.gff)
	tpms = gettpmtable(args.tpmtable)
	totaltpms = gettotaltpms(tpms)
	psis = getpsis(posfactors, tpms, totaltpms)
	with open(args.outfile, 'w') as f:
		f.write(('\t').join(['gene', 'WTSomaA', 'WTSomaB', 'WTSomaC', 'WTNeuriteA', 'WTNeuriteB', 'WTNeuriteC', 'KOSomaA', 'KOSomaB', 'KOSomaC', 'KONeuriteA', 'KONeuriteB', 'KONeuriteC']) + '\n')
		for gene in psis:
			if len(psis[gene]) != 12:
				continue
			f.write(('\t').join([gene, str(psis[gene]['WTSomaA']), str(psis[gene]['WTSomaB']), str(psis[gene]['WTSomaC']), 
				str(psis[gene]['WTNeuriteA']), str(psis[gene]['WTNeuriteB']), str(psis[gene]['WTNeuriteC']), 
				str(psis[gene]['KOSomaA']), str(psis[gene]['KOSomaB']), str(psis[gene]['KOSomaC']),
				str(psis[gene]['KONeuriteA']), str(psis[gene]['KONeuriteB']), str(psis[gene]['KONeuriteC'])]) + '\n')

