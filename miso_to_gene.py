#!/usr/bin/python

import sys
import itertools
import os
from collections import defaultdict
import tabix
from scipy.stats import fisher_exact
import numpy as np

def ensp2ensg():
	convert = defaultdict(bool)
	for line in open('/net/afterthefact/data/jmerkin/mm9_ensg.ensp.txt'):
		line = line.strip().split('\t')
		if len(line) < 2: continue
		convert[line[1]] = line[0]
	return convert

def get_go_bits_fromrob():
	go = defaultdict(set)
	all_gos = defaultdict()
	fo = open('/net/sugarman/scratch/jmerkin/annotations/mm9.genes.gff')
	fo.readline()
	for line in fo:
		line = line.strip().split('\t')
		print line[-1], line
		continue
		if line[2] != 'gene': continue
		attrs = defaultdict(bool)
		for bit in line[-1].split(';'):
			bit = bit.split('=')
			attrs[bit[0]] = bit[1]
		print attrs
		gene = attrs['ID']
		gos = attrs['GO Term Accession (bp)']
		if not (gene and gos): continue
		if not gene: continue
		for gname in gos:
			go[gene].add(gname)
			all_gos[gname] = True
	print len(all_gos)
	sys.exit()
	return go

def get_go_bits_from_david():
	go = defaultdict(set)
	names = defaultdict(bool)
	for ff in ['ENSEMBL_GENE_ID2GOTERM_BP_FAT.txt', 'ENSEMBL_GENE_ID2GOTERM_CC_FAT.txt', 
			'ENSEMBL_GENE_ID2GOTERM_MF_FAT.txt', 'ENSEMBL_GENE_ID2KEGG_PATHWAY.txt', 
			'ENSEMBL_GENE_ID2PANTHER_PATHWAY.txt', 'ENSEMBL_GENE_ID2SP_PIR_KEYWORDS.txt', 
			'ENSEMBL_GENE_ID2UP_SEQ_FEATURE.txt']:
		fbit = ff.replace('ENSEMBL_GENE_ID2','').replace('.txt', '')
		for line in open('/net/afterthefact/data/jmerkin/DAVIDKnowledgebase/mm9/'+ff):
			gene, cat = line.strip('\n').split('\t')
			go[gene].add('%s.%s' %(cat, fbit))
	
	ff = 'ENSEMBL_GENE_ID2DAVID_GENE_NAME.txt'
	for line in open('/net/afterthefact/data/jmerkin/DAVIDKnowledgebase/mm9/'+ff):
		gene, name = line.strip('\n').split('\t')
		names[gene] = name
	return go, names


def get_go_bits():
	go = defaultdict(set)
	convert = ensp2ensg()
	ensembl = 'ENSEMBL'
	GO = 'GO'
	for line in open('/net/afterthefact/data/jmerkin/gene_association.mgi'):
		if line.startswith('!'): continue
		if ensembl not in line: continue
		gene = False
		gene = False
		for bit in line.strip().split():
			bit = bit.split(':')
			if bit[0] == GO:
				gname = ':'.join(bit)
			if bit[0] != ensembl: continue
			gene = convert[bit[1]]
		if not gene: continue
		go[gene].add(gname)
	return go

global known
known = {}
def get_ate(ref, line, tab):
	def test_gene(ref, gene, isoform):
		exons = isoform.split('@')
		for exon in exons:
			for uc in exon.split('-'):
				if uc.startswith('uc'): known[uc] = True
				if uc not in ref: continue
				gene = gene or ref[uc][1]
		return gene

	gene = None

	for test in line[0].split('@')[1:]:
		if test in ref and ref[test]:
			gene = ref[test][1]
			break

	if not gene:
		isoforms = line[9].replace("'", '').split(',')
		for isos in isoforms:
			isos = isos.split('_')
			for iso in isos:
				iso = '.'.join(iso.split('.')[:-2])
				gene = test_gene(ref, gene, iso)
				if gene: break
	return gene


def test_exon(ref, chrm, start, stop):
	if (chrm, start) in ref:
		gene = ref[(chrm, start)]
		return gene
	if (chrm, stop) in ref:
		gene = ref[(chrm, stop)]
		return gene
	return None

def get_a5(ref, line, tab):
	def test_gene(gene, isoform, terms):
		exon1, exon = isoform.split('@')
		chrm, start, stop, strand = exon.split(':')
		terms.append(start)
		terms.append(stop)
		gene = test_exon(ref, chrm, start, stop)
		chrm, start, stop, strand = exon1.split(':')
		start1, start2 = stop.split('|')
		terms.append(start1)
		terms.append(start2)
		gene = gene or test_exon(ref, chrm, start1, start2)
		return gene

	gene = None
	terms = []
	gene = test_gene(gene, line[0], terms)
	if gene: return gene

	isoforms = line[9].replace("'", '').split(',')
	for isos in isoforms:
		isos = isos.split('_')
		for iso in isos:
			gene = test_gene(gene, iso.split('.')[0], terms)
			if gene: return gene
	return None

def find_by_tabix(tab, chrm, terms, strand):
	gene = None

	terms = map(int, terms)
	if len(terms) < 2: return None
	start = min(terms)
	stop = max(terms)
	genes = defaultdict(int)
	for line in tab.query(chrm.replace('chr', ''), start, stop):
		if line[6] != strand: continue
		gene = line[-1].split('"')[1]
		beg, end = map(int, line[3:5])
		thisstart = max(start, beg)
		thisstop = min(stop, end)
		genes[gene] = thisstop - thisstart
	ngenes = len(genes)
	if ngenes == 1:
		gene = genes.keys()[0]
	elif ngenes > 1:
		maxval = max(genes.values())
		these = [ii for ii in genes if genes[ii] == maxval]
		if len(these) == 1:
			gene = these[0]
	return gene

def get_tutr(ref, line, tab):
	def test_gene(gene, isoform, terms):
		for exon in isoform.split('@'):
			chrm, start, stop, strand = exon.split(':')
			terms.append(start)
			terms.append(stop)
			gene = gene or test_exon(ref, chrm, start, stop)
		return gene

	gene = None
	terms = []
	gene = test_gene(gene, line[0], terms)
	chrm, start, stop, strand = line[0].split(':')[:4]
	strand = strand[0]

	if not gene:
		isoforms = line[9].replace("'", '').split(',')
		for isos in isoforms:
			isos = isos.split('_')
			for iso in isos:
				gene = gene or test_gene(gene, iso.split('.')[0], terms)
				if gene: break
	gene = gene or find_by_tabix(tab, chrm, terms, strand)
	return gene

def get_a3(ref, line, tab):
	def test_gene(gene, isoform, terms):
		exon, exon1 = isoform.split('@')
		chrm, start, stop, strand = exon.split(':')
		terms.append(start)
		terms.append(stop)
		gene = test_exon(ref, chrm, start, stop)
		chrm, start, stop, strand = exon1.split(':')
		start1, start2 = start.split('|')
		terms.append(start1)
		terms.append(start2)
		gene = gene or test_exon(ref, chrm, start1, start2)
		return gene

	gene = None
	terms = []
	gene = test_gene(gene, line[0], terms)
	if gene: return gene

	isoforms = line[9].replace("'", '').split(',')
	for isos in isoforms:
		isos = isos.split('_')
		for iso in isos:
			gene = test_gene(gene, iso.split('.')[0], terms)
			if gene: return gene
	return None
	terms = map(int, terms)
	start = min(terms)
	stop = max(terms)
	genes = []
	for exon in line[0].split('@'):
		chrm, start, stop, strand = exon.split(':')
		start, stop = map(int, [start, stop])
		for ex in tab.query(chrm.strip('chr'), start, stop):
			print ex
			if ex[6] == strand:
				pass

	return gene

def get_as(ref, line, tab):
	gene = None
	if "ENS" in line[0]:
		return None
	for exon in line[0].split('@'):
		chrm, start, stop, strand = exon.split(':')
		if (chrm, start) in ref:
			gene = ref[(chrm, start)]
			return gene
			break
		if (chrm, stop) in ref:
			gene = ref[(chrm, stop)]
			return gene
			break
	genes = []
	for exon in line[0].split('@'):
		chrm, start, stop, strand = exon.split(':')
		start, stop = map(int, [start, stop])
		for ex in tab.query(chrm.strip('chr'), start, stop):
			if ex[6] == strand:
				pass

	return gene

def get_reffunc(misotype, known2ens, ref):
	if misotype == 'README' or 'genes' in misotype or 'analyses' in misotype: 
		this_ref = None
		reffunc = None
	elif misotype in ('ALE', 'AFE'):
		this_ref = known2ens
		reffunc = get_ate
		#reffunc = None
	elif 'Tandem' in misotype:
		this_ref = ref
		reffunc = get_tutr
		#reffunc = None
	elif 'A3' in misotype:
		this_ref = ref
		reffunc = get_a3
		#reffunc = None
	elif 'A5' in misotype:
		this_ref = ref
		reffunc = get_a5
		#reffunc = None
	elif 'RI' in misotype:
		this_ref = ref
		reffunc = get_as
		#reffunc = None
	elif 'MXE' in misotype:
		this_ref = ref
		reffunc = get_as
		#reffunc = None
	elif 'SE' in misotype:
		this_ref = ref
		reffunc = get_as
	else:
		this_ref, reffunc = None, None
	return reffunc, this_ref


def get_ref(reffile):
	ref = defaultdict(bool)
	for line in open(reffile):
		chrm, b1, ctype, start, stop, b3, strand, b4, attr = line.strip('\n').split('\t')
		#if ctype != 'exon': continue
		chrm = chrm.strip('chr')
		chrm = 'chr%s' %(chrm)
		gene = attr.split('"')[1]
		ref[(chrm, start)] = gene
		ref[(chrm, stop)] = gene
	return ref

def get_known2ens(species='mm9'):
	known2ens = defaultdict(bool)
	for line in open('/net/afterthefact/data/jmerkin/%s_known2ensembl' %(species)):
		known, enst, ensg, symbol = line.split('\t')
		known2ens[known] = enst, ensg
	for line in open('/net/afterthefact/data/jmerkin/%s_ensembl2known' %(species)):
		ensg, enst, known = line.strip('\n').split('\t')
		if not known: continue
		known2ens[known] = enst, ensg
	return known2ens

def main():
	reffile = sys.argv[1]
	indir = sys.argv[2]
	minbf = float(sys.argv[3])
	minpsi = .25
	#gene2go = get_go_bits()
	gene2go, gene2name = get_go_bits_from_david()
	go2gene = defaultdict(set)
	these, alls = [{} for _ in xrange(2)]
	faileds, types = [[] for _ in xrange(2)]
	go_fore, go_back = [defaultdict(int) for _ in xrange(2)]
	tab = tabix.Tabix('/net/afterthefact/data/jmerkin/Mus_musculus.NCBIM37.67.gtf.gz')

	known2ens = get_known2ens()
	ref = get_ref(reffile)

	conv_file = open('convert_allevents_id2gene', 'w')
	for misotype in os.listdir(indir):
		reffunc, this_ref = get_reffunc(misotype, known2ens, ref)
		if reffunc is None: continue
		types.append(misotype)
		fi = open('%s/Comparisons/N2ASoma_vs_N2AAxon/bayes-factors/N2ASoma_vs_N2AAxon.miso_bf' %(misotype, ), 'r')
		line = fi.readline()
		conv_file.write(line)
		#import code ; code.interact(local=locals())
		for line in fi:
			line = line.split('\t')
			gene = reffunc(this_ref, line, tab)
			if gene:
				if float(line[8]) > minbf and abs(float(line[7])) > minpsi:
					these[gene] = True
				alls[gene] = True
				conv_file.write('%s\t%s\t%s' %(misotype, gene, '\t'.join(line)))
			else:
				faileds.append(misotype)
	conv_file.close()
	types = '.'.join(types)

	for gene in alls:
		for go in gene2go[gene]:
			go_back[go] += 1
	for gene in these:
		for go in gene2go[gene]:
			go_fore[go] += 1
			go2gene[go].add('%s:%s' %(gene, gene2name[gene]))

	fore = len(these)
	back = len(alls)
	scoreds, folds = [[] for _ in xrange(2)]
	for go in go_back:
		back_with = go_back[go]
		fore_with = go_fore[go]
		if min(back_with, fore_with) < 2: continue
		#fore_with = max(fore_with-1 , 0)
		back_without = back - back_with 
		fore_without = fore - fore_with
		try:
			fold = float(fore_with * back) / float(back_with * fore)
		except:
			fold = 'NA'
		table = [
			#[back_with, back_without],
			#[fore_with, fore_without]
			#[back_with, fore_with],
			#[back_without, fore_without]
			[fore_with, fore_without],
			[back_with, back_without]
				]
		#print table
		table = np.array(table)
		pval = fisher_exact(table, alternative='greater')[1]
		scoreds.append((pval, go, fold))
	scoreds.sort(key=lambda xx: xx[0])
	scores, names, folds = zip(*scoreds)
	names = np.array(names)
	scores = np.array(scores)
	nscores = scores.shape[0]
	folds = np.array(folds)
	bonferroni = np.minimum(scores * float(nscores), 1.)
	benjamini = []
	oldp, knum, store = 0, 0, 0
	for ii in scores:
		benjamini.append(ii * nscores / (nscores - knum) )
		#benjamini.append(ii * (nscores - knum) / nscores )
		store += 1
		if oldp == ii:
			# to handle ties. count number of tied scores, then add them later
			pass
		else:
			knum += store
			store = 0
		oldp = ii
	benjamini = np.minimum(np.array(benjamini), 1.)

	print faileds
	print 'failed', len(faileds)
	print 'these', len(these)
	print 'all', len(alls)

	nscores = scores.shape[0]
	iis = np.arange(nscores) + 1
	Q = 0.05
	Qs = iis * Q / nscores
	def test_ben(pv, ii, ll, Q=0.05):
		if pv < ii * Q / ll:
			return True
		else:
			return False

	f_end = 'bf%s_psi%s_%s' %(minbf, minpsi, types)
	outf = open('go_analyses_%s' %(f_end), 'w')
	outf.write('term\tp-value\tbenjamini\tfdr\tfold_enrich\tgenes\n')
	passed = False
	for go, pvalue, bonf, benj, fold, ind, qv in reversed(
		zip(names, scores, bonferroni, benjamini, folds, iis, Qs)):
		if fold < 1: continue
		#if benj > 0.05: break
		#print pvalue, qv
		if passed or pvalue < qv: #test_ben(pvalue, ind, nscores, Q=Q):
			#print pvalue, benj, bonf
			line = '\t'.join(map(str, [go, pvalue, benj, qv, fold, ';'.join(go2gene[go])]))
			outf.write(line)
			outf.write('\n')
			passed = True
	outf.close()
		


	fout = open('genes_sigdif_%s' %(f_end), 'w')
	for gene in these:
		fout.write(gene)
		fout.write('\n')
	fout.close()

	fout = open('genes_all_%s' %(f_end), 'w')
	for gene in alls:
		fout.write(gene)
		fout.write('\n')
	fout.close()
			


if __name__ == '__main__':
	main()
