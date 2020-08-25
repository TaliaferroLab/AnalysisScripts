import gffutils
import argparse
import os

def split3UTR(UTR3gff, fragsize, outfile):
	gff_fn = UTR3gff
	print 'Indexing gff...'
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	UTR3s = db.features_of_type('UTR3')

	outfh = open(outfile, 'w')

	for UTR3 in UTR3s:
		#Only going to consider single exon UTRs
		if len(list(db.children(UTR3, featuretype = 'exon', level = 1))) > 1:
			continue

		ID = UTR3.attributes['ID'][0]
		parent = UTR3.attributes['Parent'][0]
		gene_id = UTR3.attributes['gene_id'][0]

		coord = UTR3.start
		counter = 1
		while coord <= UTR3.end:
			windowstart = coord
			windowend = coord + fragsize
			idfield = 'ID=' + ID + '.utr3fragment{0}'.format(counter) + ';Parent=' + parent + ';gene_id=' + gene_id
			with open(outfile, 'a') as outfh:
				outfh.write(('\t').join([str(UTR3.chrom), 'longest3UTRfrags', 'UTR3frag', str(windowstart), str(windowend), '.', str(UTR3.strand), '.', idfield]) + '\n')
			coord = coord + fragsize + 1
			counter +=1

	os.remove(db_fn)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Input gff of 3\' UTRs to split.')
	parser.add_argument('--fragsize', type = int, help = 'Fragment size for UTR fragments.')
	parser.add_argument('--output', type = str, help = 'Output gff file.')
	args = parser.parse_args()

	split3UTR(args.gff, args.fragsize, args.output)
