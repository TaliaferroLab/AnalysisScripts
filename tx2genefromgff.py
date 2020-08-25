import gffutils
import sys
import os


def gettx2gene(gff, outfile):
	tx2gene = {} #{tx : gene}
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')

	i = 0
	for gene in genes:
		i +=1
		if i % 5000 == 0:
			print('Gene {0}...'.format(i))
		genename = str(gene.id)
		for tx in db.children(gene, featuretype = 'transcript', level = 1):
			txname = str(tx.id)
			tx2gene[txname] = genename

	print('Found {0} genes for {1} transcripts.'.format(len(set(tx2gene.values())), len(tx2gene)))

	with open(outfile, 'w') as outfh:
		outfh.write('#name' + '\t' + '#name2' + '\n')
		for tx in tx2gene:
			outfh.write(tx + '\t' + tx2gene[tx] + '\n')


gettx2gene(sys.argv[1], sys.argv[2])