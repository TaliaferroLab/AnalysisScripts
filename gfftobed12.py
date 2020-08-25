import gffutils
import os
import sys


def gfftobed12(gff, outfile):
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')

	with open(outfile, 'w') as outfh:
		for gene in genes:
			for tx in db.children(gene, featuretype = 'transcript', level = 1):
				txbed12 = db.bed12(tx, block_featuretype = 'exon', thick_featuretype = 'CDS', thin_featuretype = None, name_field = 'ID', color = None)
				outfh.write(txbed12 + '\n')

gfftobed12(sys.argv[1], sys.argv[2])