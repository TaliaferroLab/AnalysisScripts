#python3

#Given a gff of the longest 'feature' for every gene (5UTR, CDS, or 3UTR), get conservation data about that feature

import gffutils
import os
import numpy as np
import pysam
import argparse

def maketempgff(gff, gene):
	#For a given gene, we want to get the exonic coordinates and look at the conservation.
	#So to do that with bedtools intersect, we can make a temporary gff that only consists of the exonic entries for that gene.

	db_fn = os.path.abspath(gff) + '.db'
	if not os.path.isfile(db_fn):
		print('GFF db does not exist!')
		sys.exit()

	db = gffutils.FeatureDB(db_fn)
	with open('temp.gff', 'w') as outfh:
		for exon in db.children(gene, featuretype = 'exon'):
			outfh.write(str(exon) + '\n')


def intersect(gff, tbx):
	phastconsvalues = []
	coords = [] #nested list of [chrm, start, stop] for each line in tempgff
	utrexonlengths = [] #lengths of all utr exons

	with open(gff, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			chrm, start, stop = line[0], int(line[3]), int(line[4])
			coords.append([chrm, start, stop])

	for coord in coords:
		utrexonlengths.append(coord[2] - coord[1])
		for row in tbx.fetch(coord[0], coord[1], coord[2], parser = pysam.asBed()):
			score = float(row.score)
			phastconsvalues.append(score)

	medphastcons = np.median(phastconsvalues)

	#Check to see if we had scores for at least some fraction of the exonic nt
	utrlength = sum(utrexonlengths)

	if len(phastconsvalues) >= (utrlength * 0.5):
		return medphastcons
	else:
		return None

def iterategff(gff, phastconsbed):
	#Iterate through a "longest feature" gff, going through the genes one by one
	#and intersecting them with a bed of phastcons values
	scores = {} #{gene : score}
	
	tbx = pysam.TabixFile(phastconsbed) #phastconsbed has to have a tabix index in the same directory

	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')

	genecount = 0
	for gene in genes:
		genecount +=1
		if genecount % 100 == 0:
			print('Gene {0}...'.format(genecount))

		#Make a gff that is only the exonic coordinates of this gene
		#This will be saved as temp.gff
		maketempgff(gff, str(gene.id))
		#Intersect these coordinates with the phastcons coordinates
		medphastcons = intersect('temp.gff', tbx)
		os.remove('temp.gff')
		scores[str(gene.id).split(':')[1].split('.')[0]] = medphastcons

	outfile = gff + '.phyloPmedians.txt'
	with open(outfile, 'w') as outfh:
		outfh.write(('\t').join(['gene', 'medianphastconsscore']) + '\n')
		for gene in scores:
			outfh.write(('\t').join([gene, str(scores[gene])]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Longest feature gff.')
	parser.add_argument('--phastconsbed', type = str, help = 'Phastcons values in bed format. Can be compressed.')
	args = parser.parse_args()

	iterategff(args.gff, args.phastconsbed)

