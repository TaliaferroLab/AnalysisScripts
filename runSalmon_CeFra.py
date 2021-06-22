from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess
import argparse
import gzip
import sys
import os

#python2


def makeindex(transcripts, k):
	#Transcripts is fasta file of transcripts to be quantified.
	#Can be gzipped.

	k = str(k)

	command = ['salmon', 'index', '-t', transcripts, '-i', 'transcripts.idx', '--type', 'quasi', '-k', k]
	subprocess.call(command)

def runsalmon(outputname, threads, reads1, reads2):

	#paired end
	if reads2:
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', outputname, '--index', 'transcripts.idx']

	#Single end, setup here for ribosome footprint data
	elif not reads2:
		fldMean = '35' #fragment length distribution mean
		fldSD = '10' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', outputname, '--index', 'transcripts.idx']

	subprocess.call(command)

def fastqlengthfilter(allowedlengths, infile):
	#Maybe you only want to allow reads of certain lengths
	#This is currently only set up to handle single end reads.  It's meant for something like ribosome footprinting.

	allowedlengths = [int(l) for l in allowedlengths]

	counter = 0
	passingreads = 0
	outfilename = 'temp.fastq'
	with gzip.open(infile, 'rb') as infh, open(outfilename, 'w') as outfh:
		try:
			for title, seq, qual in FastqGeneralIterator(infh):
				counter +=1
				if counter % 100000000 == 0:
					print 'On read {0} of {1}.'.format(counter, infile)
				if len(seq) in allowedlengths:
					passingreads +=1
					outfh.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

		except ValueError:
			print 'Title and second title line don\'t match for read {0}.'.format(title)

	print '{0} of {1} reads ({2} %) pass length filters.'.format(passingreads, counter, round((passingreads / float(counter)) * 100, 4))

def fastqtrimmer(threeprimetrim, forreads, revreads):
	#Maybe you want to trim the reads from the 3' end before giving them to salmon.
	#Trim <threeprimetrim> nt from the 3' end of the read
	threeprimetrim = int(threeprimetrim)

	counter = 0
	foutfilename = 'tempf.fastq'
	routfilename = 'tempr.fastq'
	with gzip.open(forreads, 'rb') as forinfh, gzip.open(revreads, 'rb') as revinfh, open(foutfilename, 'w') as foroutfh, open(routfilename, 'w') as revoutfh:
		try:
			for title, seq, qual in FastqGeneralIterator(forinfh):
				counter +=1
				if counter % 1000000 == 0:
					print 'On read {0} of {1}.'.format(counter, forreads)
				foroutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))
		except ValueError:
			pass

		try:
			for title, seq, qual in FastqGeneralIterator(revinfh):
				counter +=1
				if counter % 1000000 == 0:
					print 'On read {0} of {1}.'.format(counter, forreads)
				revoutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))

		except ValueError:
			pass

	print 'Done trimming {0} and {1}.'.format(forreads, revreads)





if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--transcripts', type = str, help = 'Fasta file of transcripts to quantify.')
	parser.add_argument('-k', type = int, help = 'Value of k to use when making the index.')
	parser.add_argument('--allowedlengths', type = str, help = 'Optional. Comma separated list of allowed read lengths.')
	parser.add_argument('--threeprimetrim', type = int, help = 'Optional. Number of bases to trim from the 3\' end of reads.')
	parser.add_argument('--labels', type = str, help = 'Comma separated list of labels for experiments.')
	parser.add_argument('--threads', type = str)
	parser.add_argument('--forwardreads', type = str, help = 'Comma separated list of forward fastq files in same order as labels.')
	parser.add_argument('--reversereads', type = str, help = 'Comma separated list of reverse fastq files in same order as labels. Only needed with paired end reads.')
	args = parser.parse_args()

	print 'Indexing transcripts...'
	makeindex(args.transcripts, args.k)
	print 'Done indexing!'

	rawreaddir = '/beevol/home/taliaferro/data/CeFra/RawReads'
	for rep in ['Rep1', 'Rep2']:
		for cellline in ['HepG2', 'K562']:
			for libprep in ['polyA', 'ribodep']:
				for location in ['cytosol', 'insoluble', 'membrane', 'nucleus', 'total']:
					fr = '{0}.{1}.{2}.{3}.1.fastq.gz'.format(cellline, location, libprep, rep)
					rr = '{0}.{1}.{2}.{3}.2.fastq.gz'.format(cellline, location, libprep, rep)
					label = '{0}_{1}_{2}_{3}'.format(cellline, location, libprep, rep)
					freads = os.path.join(rawreaddir, fr)
					rreads = os.path.join(rawreaddir, rr)
					runsalmon(label, args.threads, freads, rreads)

#	for rep in ['Rep1', 'Rep2', 'Rep3', 'Rep4', 'Rep5', 'Rep6']:
#		fr = 'KDEL-HRP_Input_{0}_1.fastq.gz'.format(rep)
#		rr = 'KDEL-HRP_Input_{0}_2.fastq.gz'.format(rep)
#		label = 'KDEL-HRP_Input_{0}'.format(rep)
#		freads = os.path.join(rawreaddir, fr)
#		rreads = os.path.join(rawreaddir, rr)
#		runsalmon(label, args.threads, freads, rreads)

#	for samp in ['Mito', 'NES', 'NLS']:
#		for frac in ['Input', 'RIP']:
#			for rep in ['Rep1', 'Rep2', 'Rep3', 'Rep4', 'Rep5', 'Rep6']:
#				label = samp + '_' + 'APEX2' + '_' + frac + '_' + rep
#				fr = '{0}-APEX2_{1}_{2}_1.fastq.gz'.format(samp, frac, rep)
#				rr = '{0}-APEX2_{1}_{2}_2.fastq.gz'.format(samp, frac, rep)
#				freads = os.path.join(rawreaddir, fr)
#				rreads = os.path.join(rawreaddir, rr)
#				runsalmon(label, args.threads, freads, rreads)
