#Given a gff annotation, make a fasta for all transcript sequences.

import gffutils
import os
from Bio import SeqIO
import argparse
import gzip

#tags we dont want
#cds_end_NF, cds_start_NF, mRNA_end_NF, mRNA_start_NF, TSL = 4, TSL = 5

def filtertranscripts(transcript, db):
	attributepass = False
	tslpass = False

	atts = transcript.attributes
	if 'tag' not in atts or ('mRNA_end_NF' not in atts['tag'] and 'mRNA_start_NF' not in atts['tag'] and 'cds_start_NF' not in atts['tag'] and 'cds_start_NF' not in atts['tag']):
		attributepass = True

	try:
		tsl = transcript.attributes['transcript_support_level'][0]
		if tsl == 'NA':
			pass
		elif int(tsl) <= 3:
			tslpass = True
	except KeyError:
		pass

	if attributepass and tslpass:
		return True
	else:
		return False



def getcDNAcoords(gff, usefilters):
	cDNAcoords = {} #{ENSTRANS_chrm_strand: [[exon1start, exon1stop], [exon2start, exon2stop]]}
	genecount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1 
		if genecount % 10000 == 0:
			print 'Gene {0}...'.format(genecount)

		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start', level = 1):
			transcriptpassesfilters = filtertranscripts(transcript, db)
			if usefilters and not transcriptpassesfilters:
				continue

			transcriptname = str(transcript.id) + '_' + str(transcript.chrom) + '_' + str(transcript.strand)
			exoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
			for exon in db.children(transcript, featuretype = 'exon', order_by = 'start', level = 1):
				exoncoords.append([exon.start, exon.end])
			if exoncoords:
				cDNAcoords[transcriptname] = exoncoords

	return cDNAcoords

def cDNAcoordstoseq(cDNAcoords, genomefasta, outputfasta):
	cDNAseqs = {} #{ENSTRANS : seq}
	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'
	chrmswithoutseq = [] #Chromsome names that are in allCDScoords but that don't have a fasta entry in genomefasta

	for transcript in cDNAcoords:
		cDNAseq = ''
		#Ensembl vs UCSC transcript names
		#Ensembl
		if 'NM_' not in transcript and 'NR_' not in transcript:
			tname = transcript.split('_')[0].replace('transcript:', '') #in cDNAcoords txnames are like 'transcript:ENSMUST00000038375_chr2_+'
			chrm = transcript.split('_')[1]
			strand = transcript.split('_')[2]
		
		#UCSC
		elif 'NM_' in transcript or 'NR_' in transcript:
			if ':' in transcript:
				tname = transcript.split('_chr')[0].split(':')[1] #transcript:NM_175684_chr18_-
			else:
				tname = transcript.split('_chr')[0]
			chrm = transcript.split('_')[-2]
			strand = transcript.split('_')[-1]


		#Is this chromosome in genomefasta?
		if chrm not in seq_dict:
			if chrm not in chrmswithoutseq:
				print 'WARNING: No entry for chromosome {0} in genomefasta.'.format(chrm)
				print transcript
				chrmswithoutseq.append(chrm)
			continue

		for exon in cDNAcoords[transcript]:
			start = exon[0]
			end = exon[1]
			if strand == '+':
				exonseq = seq_dict[chrm].seq[start-1:end].upper()
				cDNAseq += exonseq
			elif strand == '-':
				exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
				newseq = exonseq + cDNAseq
				cDNAseq = newseq

		cDNAseqs[tname] = cDNAseq

	print 'Found sequence data for {0} of {1} transcripts.'.format(len(cDNAseqs), len(cDNAcoords))

	with open(outputfasta, 'w') as f:
		for tname in cDNAseqs:
			seq = str(cDNAseqs[tname])

			#If transcript IDs have a dot (like ENSMUST0000046543.2)
			if '.' in tname:
				tname = tname.split('.')[0]
			f.write('>' + tname + '\n' + seq + '\n')

	return cDNAseqs


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF annotation of transcripts.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--outputfasta', type = str, help = 'Output file for cDNA sequences.')
	parser.add_argument('--usefilters', action = 'store_true', help = 'Do you want all transcripts or only those that have reasonable evidence?')
	args = parser.parse_args()

	cDNAcoords = getcDNAcoords(args.gff, args.usefilters)
	cDNAcoordstoseq(cDNAcoords, args.genomefasta, args.outputfasta)


