from Bio import SeqIO
import gffutils
import gzip

def getdistinctregions(gff, genomefasta):
	distinctregions = {} #{geneid : {transcriptid(s) : [3UTR number, distinctUTRseq]}}
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		distinctseqs = {} #{transcriptid(s) : [pAsite counter (may be different than number of UTRs because not all UTRs are represented here, distinctUTRseq]}
		seenseqs = []
		utrcounter = 0
		mostdownstreamcoord = 0 #The most downstream coordinate of any UTR we've seen so far for this gene.
		geneid = str(gene.id).replace('gene:', '').split('.')[0]
		if gene.strand == '+':
			for UTR3 in db.children(gene, featuretype = 'UTR3', level = 1, order_by = 'end'):
				distinctseq = ''
				UTRid = str(UTR3.id).replace('UTR3:', '').split('.')[0]

				#If this is the first UTR for this gene
				if utrcounter == 0:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'start'):
						exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].upper()
						distinctseq += exonseq
					mostdownstreamcoord = UTR3.end
					utrcounter +=1
					distinctseqs[UTRid] = [utrcounter, str(distinctseq)]
				elif utrcounter >= 1:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'start'):
						#If this exon is somehow contained within the last one (should not be possible), skip it
						if exon.end <= mostdownstreamcoord:
							pass
						elif exon.end > mostdownstreamcoord:
							if exon.start < mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[mostdownstreamcoord:exon.end].upper()
								distinctseq += exonseq
							elif exon.start >= mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start - 1:exon.end].upper()
								distinctseq += exonseq
					
					mostdownstreamcoord = UTR3.end

					#Only going to call something a new polyA site if it's at least 50 nt away from the previous one
					#As a proxy for this, it must have at least 50 nt of "distinct" sequence
					if len(str(distinctseq)) >= 50:
						utrcounter +=1
						distinctseqs[UTRid] = [utrcounter, str(distinctseq)]

		elif gene.strand == '-':
			for UTR3 in db.children(gene, featuretype = 'UTR3', level = 1, order_by = 'start', reverse = True):
				distinctseq = ''
				UTRid = str(UTR3.id).replace('UTR3:', '').split('.')[0]

				#If this is the first UTR for this gene
				if utrcounter == 0:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'end', reverse = True):
						exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].reverse_complement().upper()
						#Must prepend instead of append this time
						distinctseq = distinctseq + exonseq
					mostdownstreamcoord = UTR3.start
					utrcounter +=1
					distinctseqs[UTRid] = [utrcounter, str(distinctseq)]
				elif utrcounter >= 1:
					for exon in db.children(UTR3, featuretype = 'exon', level = 1, order_by = 'end', reverse = True):
						#If this exon is somehow contained within the last one (should not be possible), skip it
						if exon.start >= mostdownstreamcoord:
							continue
						elif exon.start < mostdownstreamcoord:
							if exon.end > mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start-1:mostdownstreamcoord-1].reverse_complement().upper()
								distinctseq = distinctseq + exonseq
							elif exon.start <= mostdownstreamcoord:
								exonseq = seq_dict[exon.chrom].seq[exon.start-1:exon.end].reverse_complement().upper()
								distinctseq = distinctseq + exonseq


					mostdownstreamcoord = UTR3.start
					if len(str(distinctseq)) >= 50:
						utrcounter +=1
						distinctseqs[UTRid] = [utrcounter, str(distinctseq)]

		distinctregions[geneid] = distinctseqs

	return distinctregions