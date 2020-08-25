import gffutils
import os
import sys
import subprocess

#python3

#The point of this is to make a gff for use with QAPA that is as similar to the genes/transcripts that LABRAT sees with all of its filters.
#So we are going to use the same filters for this that we normally use with LABRAT
#This will take a gff as an input and output a gtf

#Does a gene pass filters?
def genefilters(gene, db):
	proteincoding = True #don't want this filter

	if 'protein_coding' in gene.attributes['gene_type']:
		proteincoding = True

	if proteincoding == True:
		return True
	else:
		return False


#Does a transcript pass filters?
def transcriptfilters(transcript, db):
	exonnumberpass = False
	TFlengthpass = False
	proteincoding = False #don't want this filter #turned on for now
	mrnaendpass = False
	#How many exons does it have
	if len(list(db.children(transcript, featuretype = 'exon'))) >= 2:
		exonnumberpass = True
	else:
		return False
	
	#What is the length of the terminal fragment
	exons = []
	if transcript.strand == '+':
		for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
			exons.append([exon.start, exon.end + 1])
	elif transcript.strand == '-':
		for exon in db.children(transcript, featuretype = 'exon', order_by = 'start', reverse = True):
			exons.append([exon.start, exon.end + 1])
	penultimateexonlen = len(range(exons[-2][0], exons[-2][1]))
	lastexonlen = len(range(exons[-1][0], exons[-1][1]))
	TFlength = penultimateexonlen + lastexonlen
	if TFlength > 200:
		TFlengthpass = True

	#Is this transcript protein coding
	if 'protein_coding' in transcript.attributes['transcript_type']:
		proteincoding = True

	#Are we confident in the 3' end of this mrnaendpass
	if 'tag' not in transcript.attributes or 'mRNA_end_NF' not in transcript.attributes['tag']:
		mrnaendpass = True

	if exonnumberpass and TFlengthpass and proteincoding and mrnaendpass:
		return True
	else:
		return False

def filtergff(gff, filteredoutfile):

	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	with open(filteredoutfile, 'w') as outfh:
		genes = db.features_of_type('gene')
		for gene in genes:
			passgenefilters = genefilters(gene, db)
			if passgenefilters == False:
				continue

			passingtranscriptids = []
			for transcript in db.children(gene, featuretype = 'transcript', level = 1):
				passtranscriptfilters = transcriptfilters(transcript, db)
				if passtranscriptfilters == True:
					passingtranscriptids.append(transcript.id)

			#If this gene has passing transcripts, write the gene
			if passingtranscriptids:
				outfh.write(str(gene) + '\n')
				for transcriptid in passingtranscriptids:
					transcript = db[transcriptid]
					outfh.write(str(transcript) + '\n')
					for exon in db.children(transcript, featuretype = 'exon', level = 1):
						outfh.write(str(exon) + '\n')
					for CDS in db.children(transcript, featuretype = 'CDS', level = 1):
						outfh.write(str(CDS) + '\n')
					for startcodon in db.children(transcript, featuretype = 'start_codon', level = 1):
						outfh.write(str(startcodon) + '\n')
					for stopcodon in db.children(transcript, featuretype = 'stop_codon', level = 1):
						outfh.write(str(stopcodon) + '\n')
					for UTR in db.children(transcript, featuretype = 'UTR', level = 1):
						outfh.write(str(UTR) + '\n')

#Ok now that we have a gff that is filtered like the one that LABRAT uses, convert that to a gff using gffread
def rungffread(ingff, outgtf):
	command = ['gffread', '-E', '-O', '-T', ingff, '-o', outgtf]
	subprocess.call(command)



filtergff(sys.argv[1], 'filtered.gff')
rungffread('filtered.gff', sys.argv[2])
os.remove('filtered.gff')
