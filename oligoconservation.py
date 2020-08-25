import pysam
import numpy as np
import argparse
import gffutils
import sys
import os


def getoligoscore(tbx, chrm, start, stop):
	scores = []

	for row in tbx.fetch(chrm, start, stop, parser = pysam.asBed()):
		score = float(row.score)
		scores.append(score)

	return len(scores), scores

def iterateoligos(gff, scorebed, outfile):
	tbx = pysam.TabixFile(scorebed) #scorebed has to have a tabix index in the same directory

	scoresdict = {} #{oligoid : [[chrm, start, stop], number of nts that have a score, median score]}
	oligocounter = 0
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	oligos = db.features_of_type('oligo')

	oligocounter = 0
	for oligo in oligos:
		oligocounter +=1
		if oligocounter % 1000 == 0:
			print('Analyzing oligo {0}...'.format(oligocounter))

		oligoid = str(oligo.id)
		genename = oligo.attributes['gene_name'][0]
		chrm = str(oligo.chrom)
		if oligo.attributes['oligo_type'][0] == 'regular_oneexon':
			start = oligo.start
			end = oligo.end
			#If this is the first oligo of a UTR, it contains 44 nt of coding sequence.  Remove it.
			if oligoid[-2:] == '.1':
				if oligo.strand == '+':
					start = oligo.start + 44
				elif oligo.strand == '-':
					end = oligo.end - 44

			ntwithscore, scores = getoligoscore(tbx, chrm, start, end)
			oligoscores = scores

		elif oligo.attributes['oligo_type'][0] == 'junction':
			junction_ntwithscore = 0
			jpscores = [] #[[scores of jp1], [scores of jp2], ...]
			for junctionpiece in db.children(oligo, featuretype = 'junctionpiece'):
				jpstart = junctionpiece.start
				jpend = junctionpiece.end
				if oligoid[-2:] == '.1' and junctionpiece.id[-2:] == '.1':
					if oligo.strand == '+':
						if jpstart + 44 < jpend:
							jpstart = jpstart + 44
						else:
							continue
					elif oligo.strand == '-':
						if jpend - 44 > jpstart:
							jpend = jpend - 44
						else:
							continue
				try:
					ntwithscore, scores = getoligoscore(tbx, chrm, jpstart, jpend)
				except ValueError: #junctionpiece was less than 44 nt long
					print(oligoid)
					continue
				junction_ntwithscore += ntwithscore
				jpscores += scores
			oligoscores = jpscores

		elif oligo.attributes['oligo_type'][0] == 'regular_multiexon':
			oligo_ntwithscore = 0
			opscores = [] #[[scores of op1], [scores of op2], ...]
			for oligopiece in db.children(oligo, featuretype = 'oligopiece'):
				opstart = oligopiece.start
				opend = oligopiece.end
				if oligoid[-2:] == '.1' and oligopiece.id[-2:] == '.1':
					if oligo.strand == '+':
						if opstart + 44 < opend:
							opstart = opstart + 44
						else:
							continue
					elif oligo.strand == '-':
						if opend - 44 > opstart:
							opend = opend - 44
						else:
							continue
				try:		
					ntwithscore, scores = getoligoscore(tbx, chrm, opstart, opend)
				except ValueError: #oligopiece was less than 44 nt long
					print(oligoid)
					continue
				oligo_ntwithscore += ntwithscore
				opscores += scores
			oligoscores = opscores

		else:
			print(oligo.attributes['oligo_type'])
			sys.exit()

		scoresdict[oligoid + '|' + genename] = [[chrm, str(oligo.start), str(oligo.stop)], str(ntwithscore), str(np.median(oligoscores)), oligoscores]

	with open(outfile, 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'chrm', 'start', 'stop', 'ntwithscore', 'medianscore', 'ntscore']) + '\n')
		for oligo in scoresdict:
			oligoscores = scoresdict[oligo][3]
			for ntscore in oligoscores:
				chrm = scoresdict[oligo][0][0]
				start = scoresdict[oligo][0][1]
				stop = scoresdict[oligo][0][2]
				ntwithscore = scoresdict[oligo][1]
				medianscore = scoresdict[oligo][2]
				outfh.write(('\t').join([oligo, chrm, start, stop, ntwithscore, medianscore, str(ntscore)]) + '\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF of oligos.')
	parser.add_argument('--scorebed', type = str, help = 'Bed file of conservation scores.  Must have tabix index in same directory.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	iterateoligos(args.gff, args.scorebed, args.outfile)