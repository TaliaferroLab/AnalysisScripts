#Get frame enrichments for ribosome footprinting data.
#Bams are 0-based. Sam is 1-based.

import gffutils
import os
import sys
import pysam
import subprocess
import argparse

#Filter bam for only 30 or 31 nt reads
def filterbam(bam, lengthfilter):
	readcounter = 0
	filteredreadcounter = 0
	with pysam.AlignmentFile(bam, 'rb') as inf:
		header = inf.header

		with pysam.AlignmentFile('temp.bam', 'wb', header = header) as outf:
			for read in inf.fetch(until_eof = True):
				readcounter +=1
				if read.query_length == lengthfilter:
					filteredreadcounter +=1
					outf.write(read)

	print 'Filtered {0} reads to find {1} ({2}%) that were {3} nt long.'.format(readcounter, filteredreadcounter, round((filteredreadcounter / float(readcounter)), 2) * 100, lengthfilter)
	return filteredreadcounter, readcounter

#Starting with a LongestCDS gff, make gff of only CDS exons.
#Then get coverage of CDS exons. Positions are 1-based.
#Output file is coverage.txt
#Not actually used here.
def getCDScoverage(bam, CDSgff):
	coveragedict = {} #{chrm : {strand : {pos : coverage}}}
	with open(CDSgff, 'r') as infh, open('temp.gff', 'w') as outfh:
		for line in infh:
			line = line.strip().split('\t')
			if line[2] == 'exon':
				outfh.write(('\t').join(line) + '\n')

	bedtoolscommand = 'bedtools coverage -a temp.gff -b {0} -d -s'.format(bam).split()
	print 'Calculating per nucleotide CDS coverage for {0}...'.format(bam)
	with open('coverage.txt', 'w') as outfh:
		subprocess.call(bedtoolscommand, stdout = outfh)

	os.remove('temp.gff')



#Given an annotation as gff, get the longest CDS for every gene.
#Then get lists of nucleotide positions for every CDS for the purpose of assigning a frame to every nt in a "longest CDS".
def getCDScoords(gff):
	allCDScoords = {} #{ENSGENE_chrm_strand : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
	allCDSnt = {} #{chrm : {strand : [[list of nucleotides of CDS1], [list of nucleotides of CDS2]]}}
	framedict = {} #{chrm : {strand : {position (1-based) : frame}}}
	geneboundaries = {} # {ensid : [genestart, genestop]}
	genecount = 0
	geneswithcodingtranscript = 0
	e2sdict = {} #{ENSGene : shortname}

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	for gene in genes:
		genecount +=1
		if genecount % 10000 == 0:
			print 'Gene {0}...'.format(genecount)
		geneboundaries[str(gene.id).replace('gene:', '')] = [gene.start, gene.end]

		#If this gene doesn't have at least one coding 'transcript' (i.e. not 'NMD_transcript_variant', and not a ncRNA), skip it
		codingtranscript = False
		for transcript in db.children(gene, featuretype = 'transcript', level = 1):
			if len(list(db.children(transcript, featuretype = 'CDS', level = 1))) > 0:
				codingtranscript = True
				geneswithcodingtranscript +=1
				break
		if not codingtranscript:
			continue


		CDSlengths = {} #{transcriptID : combined_length_of_coding_exons}
		CDScoords = {} #{transcriptID : [[cdsexon1start, cdsexon1stop], [cdsexon2start, cdsexon2stop]]}
		genename = str(gene.id)
		chrm = str(gene.chrom)
		strand = gene.strand
		for transcript in db.children(gene, featuretype = 'transcript', order_by = 'start'):
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

	#os.remove(db_fn)
	#Now reorganize allCDScoords into allCDSnt
	for gene in allCDScoords:
		chrm = gene.split('_')[1]
		strand = gene.split('_')[2]
		if chrm not in allCDSnt:
			allCDSnt[chrm] = {}
		if strand not in allCDSnt[chrm]:
			allCDSnt[chrm][strand] = []

		cdsnt = []
		for exon in allCDScoords[gene]:
			exonnt = range(exon[0], exon[1] + 1)
			cdsnt += exonnt

		#Flip the coords for minus strand
		if strand == '-':
			cdsnt = cdsnt[::-1]

		#Some CDS coords (about 2%) do not have a multiple of 3 length. They tend to be for not well-annotated genes (e.g. Gm20946)
		if len(cdsnt) %3 == 0:
			allCDSnt[chrm][strand].append(cdsnt)

	print 'Looked through {0} genes. {1} of them had non-NMD coding transcripts. Found longest CDS sequences for {2} of them.'.format(genecount, geneswithcodingtranscript, len(allCDScoords))

	cdsnt = 0
	nonframe = 0
	for chrm in allCDSnt:
		for strand in allCDSnt[chrm]:
			for cds in allCDSnt[chrm][strand]:
				if len(cds) % 3 != 0:
					nonframe +=1
				for nt in cds:
					cdsnt +=1

	#Now populate framedict
	for chrm in allCDSnt:
		if chrm not in framedict:
			framedict[chrm] = {}
		for strand in allCDSnt[chrm]:
			if strand not in framedict[chrm]:
				framedict[chrm][strand] = {}
			for CDS in allCDSnt[chrm][strand]:
				for ind, position in enumerate(CDS):
					frame = ind % 3
					framedict[chrm][strand][position] = frame

	return framedict

#Given the frame of all CDS positions and the bam, get the starting position of all reads
#and then get their frames
def getframes(bam, positionframedict):
	frames = {} #{0 : <number of reads in plus1 frame>, 1: <number of reads in plus 2 frame>, etc.}
	unfindablepositions = 0 #number of positions in coveragefile that are not present in positionframedict
	coveragefileline = 0
	for frame in [0, 1, 2]:
		frames[frame] = 0
	with pysam.AlignmentFile(bam, 'rb') as infh:
		for read in infh.fetch(until_eof = True):

			#For non-RFP reads, the first read of mate pair is on the wrong strand
			#Also, only consider unique mappers
			if read.is_read1 or read.get_tag('NH') > 1 or read.is_secondary:
				continue

			if read.is_reverse:
				strand = '-'
			elif not read.is_reverse:
				strand = '+'
			chrm = read.reference_name
			if strand == '+':
				readstart = read.reference_start + 1
			elif strand == '-':
				readstart = read.reference_end

			#Offset readstart by 12 to get P site (or 13 for 31 or 32 nt reads)
			refpositions = read.get_reference_positions(full_length = True)
			#Using full_length = True will allow soft-clipped nts to have placeholders of None in this list
			#We still want them as placeholders so that the offset will be correct
			refpositions = [pos + 1 if type(pos) == int else pos for pos in refpositions] #make 1 based
			if read.query_length == 30 or read.query_length == 31 or read.query_length == 32 or read.query_length == 33 or read.query_length == 34 or read.query_length == 35:
				if strand == '+':
					refpos = refpositions[13]
				elif strand == '-':
					refpositions.reverse()
					refpos = refpositions[13]
			else:
				if strand == '+':
					refpos = refpositions[12]
				elif strand == '-':
					refpositions.reverse()
					refpos = refpositions[12]
			#Not all readstarts are going to be in positionframedict.
			#Not all reads map to CDS.
			try:
				frame = positionframedict[chrm][strand][refpos]
				frames[frame] +=1
			except KeyError:
				pass

	totalCDSreadcount = sum(frames.values())
	#If there were no reads mapping to a CDS
	if totalCDSreadcount == 0:
		framesnormalized = {0:0, 1:0, 2:0}
		
	else:
		framesnormalized = {}
		for frame in [0,1,2]:
			normalizedcount = frames[frame] / float(totalCDSreadcount)
			framesnormalized[frame] = normalizedcount

	return framesnormalized

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--bamdirectory', type = str, help = 'Directory containing all bam files to process.')
	parser.add_argument('--gff', type = str, help = 'Gff of genome annotation.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	args = parser.parse_args()

	#readlengths = [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]
	readlengths = range(15, 51)


	positionframedict = getCDScoords(args.gff)
	with open(args.outfile, 'w') as outfh:
		outfh.write(('\t').join(['sample', 'readlength', 'lengthfraction', 'frame', 'framefraction']) + '\n')
		for f in os.listdir(args.bamdirectory):
			if f.endswith('.bam'):
				samplename = f.split('.')[0] + '.' + f.split('.')[1]
				for readlength in readlengths:
					filteredreadcount, readcount = filterbam(f, readlength)
					fraction = filteredreadcount / float(readcount)
					print 'Getting frames for {0} nt reads in {1}...'.format(readlength, f)
					frames = getframes('temp.bam', positionframedict)
					for frame in sorted(frames):
						outfh.write(('\t').join([samplename, str(readlength), str(fraction), str(frame), str(frames[frame])]) + '\n')

	os.remove('temp.bam')


