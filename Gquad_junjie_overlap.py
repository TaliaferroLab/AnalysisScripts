import gffutils
import os
import sys

#Given a gff of junjie's gquad regions and a gff of gene regions, caluclate the number of basepair overlaps for each gene region

def getjunjieregions(junjiegff):
	#Make gff databases
	print 'Indexing gff...'
	gff_fn = junjiegff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	junjieregions = {} #{chrm : {strand : [list of nt in a junjie region]}}

	regions = db.features_of_type('region')

	for region in regions:
		for exon in db.children(region, featuretype = 'exon'):
			chrm = exon.chrom
			strand = exon.strand
			nt = range(exon.start, exon.end + 1)

			if chrm not in junjieregions:
				junjieregions[chrm] = {}
			if strand not in junjieregions[chrm]:
				junjieregions[chrm][strand] = []

			junjieregions[chrm][strand] += nt

	#Remove duplicates
	for chrm in junjieregions:
		for strand in junjieregions[chrm]:
			nt = junjieregions[chrm][strand]
			junjieregions[chrm][strand] = list(set(nt))

	os.remove(db_fn)
	return junjieregions

def overlapwithtxregions(junjieregions, gff, outfile):

	#Make gff databases
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	overlaps = {} #{geneid : number of nt that overlap}

	genes = db.features_of_type('gene')

	for gene in genes:
		genent = []
		geneid = gene.attributes['gene_id'][0].split('.')[0]
		for exon in db.children(gene, featuretype = 'exon'):
			chrm = exon.chrom
			strand = exon.strand
			exonnt = range(exon.start, exon.end + 1)
			genent += exonnt

		try: #not all chrms (like chrY) are in junjieregions
			g4nt = junjieregions[chrm][strand]
			overlap = list(set(g4nt).intersection(genent))
		except KeyError:
			overlap = []
	
		overlaps[geneid] = len(overlap)

	os.remove(db_fn)

	with open(outfile, 'w') as outfh:
		outfh.write(('\t').join(['Gene', 'junjieg4overlap']) + '\n')
		for gene in overlaps:
			outfh.write(('\t').join([gene, str(overlaps[gene])]) + '\n')




regions = getjunjieregions(sys.argv[1])
overlaps = overlapwithtxregions(regions, sys.argv[2], sys.argv[3])