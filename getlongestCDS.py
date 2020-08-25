import gffutils
import os
from Bio import SeqIO
import argparse

def getCDScoords(gff):
	allCDScoords = {} #{ENSGENE_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	genecount = 0
	noCDSexonscount = 0

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge')

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1

		#If this gene has 0 coding exons (for example a ncRNA), skip it
		CDS_count = len(list(db.children(gene, featuretype = 'mRNA')))
		if CDS_count == 0:
			noCDSexonscount +=1
			continue

		CDSlengths = {} #{transcriptID : combined_length_of_coding_exons}
		CDScoords = {} #{transcriptID : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
		genename = str(gene.id)
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'mRNA', order_by = 'start'):
			transcriptID = str(transcript.id)
			CDScoords[transcriptID] = []
			CDSlength = 0
			for codingexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
				CDScoords[transcriptID].append([codingexon.start, codingexon.end])
				exonlength = codingexon.end - codingexon.start
				CDSlength += exonlength
			CDSlengths[transcriptID] = CDSlength

		longestcds = max(CDSlengths.iterkeys(), key = (lambda key: CDSlengths[key]))
		for transcript in CDScoords:
			if transcript == longestcds:
				allCDScoords[genename + '_' + chrm + '_' + strand] = CDScoords[transcript]

	os.remove(db_fn)

	print 'Looked through {0} genes. {1} of them had no coding exons. Found longest CDS sequences for {2} of them.'.format(genecount, noCDSexonscount, len(allCDScoords))
	return allCDScoords

def getCDSSequences(allCDScoords, genomefasta, ens2short):
	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	print 'Done indexing!'
	CDSseqs = {} #{genename : CDSseq}
	e2sdict = {} #{ENSGene : shortname}

	infh = open(ens2short, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0].startswith('ENSMUSG'):
			e2sdict[line[0]] = line[1]
	infh.close()

	for gene in allCDScoords:
		CDSseq = ''
		genename = gene.split('_')[0]
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]
		for coords in allCDScoords[gene]:
			start = coords[0]
			end = coords[1]
			if strand == '+':
				exonseq = seq_dict[chrm].seq[start-1:end].upper()
				CDSseq += exonseq
			elif strand == '-':
				exonseq = seq_dict[chrm].seq[start-1:end].reverse_complement().upper()
				newseq = exonseq + CDSseq
				CDSseq = newseq

		if genename in e2sdict:
			shortname = e2sdict[genename]
			CDSseqs[shortname] = str(CDSseq)

	return CDSseqs

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Gff file containing CDS coordinates.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
	parser.add_argument('--ens2short', type = str, help = 'File containing conversions between ENSG names and gene short names. Ensembl2GeneName_Mouse.txt for example.')
	parser.add_argument('--output', type = str, help = 'Output file.')
	args = parser.parse_args()




allCDScoords = getCDScoords(args.gff)
seqs = getCDSSequences(allCDScoords, args.genomefasta, args.ens2short)
outfh = open(args.output, 'w')
for gene in seqs:
	outfh.write('>' + gene + '\n')
	outfh.write(seqs[gene] + '\n')
outfh.close()
			