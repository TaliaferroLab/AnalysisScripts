#!/usr/bin/python

import os
import sys
import tabix
from collections import defaultdict
from miso_to_gene import get_go_bits_from_david, get_known2ens, get_reffunc

def get_ref(reffiles):
	ref = defaultdict(set)
	for reffile in reffiles:
		for line in open(reffile):
			chrm, b1, ctype, start, stop, b3, strand, b4, attr = line.strip('\n').split('\t')
			#if ctype != 'exon': continue
			chrm = chrm.strip('chr')
			chrm = 'chr%s' %(chrm)
			gene = attr.split('"')[1]
			ref[(chrm, start)].add(gene)
			ref[(chrm, stop)].add(gene)
	return ref

if __name__ == '__main__':

	reffiles = sys.argv[1].split(',')
	indir = sys.argv[2]

	gene2go, gene2name = get_go_bits_from_david()
	known2ens = get_known2ens()
	ref = get_ref(reffiles)
	tab = tabix.Tabix('/net/afterthefact/data/jmerkin/Mus_musculus.NCBIM37.67.gtf.gz')
	types = []
	faileds = []

	conv_file = open('convert_allevents_id2gene', 'w')

	for misotype in os.listdir(indir):
		reffunc, this_ref = get_reffunc(misotype, known2ens, ref)
		print misotype, reffunc
		if reffunc is None: continue
		types.append(misotype)
		fi = open('%s/%s/Comparisons/N2ASoma_vs_N2AAxon/bayes-factors/N2ASoma_vs_N2AAxon.miso_bf' %(indir, misotype, ), 'r')
		conv_file.write(fi.readline())
		#import code ; code.interact(local=locals())
		for line in fi:
			line = line.split('\t')
			gene = reffunc(this_ref, line, tab)
			if gene:
				gene = ','.join(gene)
				conv_file.write('%s\t%s\t%s' %(misotype, gene, '\t'.join(line)))
			else:
				conv_file.write('%s\t%s\t%s' %(misotype, None, '\t'.join(line)))
				faileds.append(misotype)
	conv_file.close()
