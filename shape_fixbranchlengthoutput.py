#Because some oligos overlap, there are some kmers whose coordinates are contained within more than one oligo.
#This is designed to take the output from Jason's branch length script and put it into a more usable format.

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def getRvalues(rvalues):
	rvaluedict = {} # {oligoname : Rvalue}
	rfh = open(rvalues, 'r')
	for line in rfh:
		line = line.strip().split('\t')
		if line[0] == 'species':
			continue
		rvaluedict[line[1]] = line[2]

	return rvaluedict
	rfh.close()

def getOligocoords(gff):
	oligodict = {} # {oligoname : [chrm, start, stop, strand, age, location]}
	gfffh = open(gff, 'r')
	for line in gfffh:
		line = line.strip().split('\t')
		chrm = line[0]
		age = line[1]
		location = line[2]
		start = int(line[3])
		stop = int(line[4])
		strand = line[6]
		oligoname = line[8][3:]
		oligodict[oligoname] = [chrm, start, stop, strand, age, location]

	print 'Read {0} oligos.\n'.format(len(oligodict))

	gfffh.close()

	return oligodict

def addfieldstobranchlength(blfile, oligodict, rvaluedict, outfile):
	blfh = open(blfile, 'r')
	outfh = open(outfile, 'w')
	outfh.write(('\t').join(['Oligoname', 'chrm','oligostart','oligostop','strand','age','location','rvalue','kmerstart','kmerstop','branchlength','kmer']) + '\n')
	bedcounter = 0
	for line in blfh:
		match = False
		bedcounter +=1
		if bedcounter % 1000 == 0:
			print 'Considering bed line {0}...'.format(bedcounter)
		line = line.strip().split('\t')
		if len(line) < 6: #some lines in the bed file seem to be missing the kmer sequence
			continue
		chrm = line[0]
		kmerstart = int(line[1]) + 1 #convert bed coords to gff coords
		kmerstop = int(line[2])
		branchlength = float(line[4])
		kmer = Seq(line[5], generic_dna)
		for oligo in oligodict:
			oligochrm = oligodict[oligo][0]
			oligostart = int(oligodict[oligo][1])
			oligostop = int(oligodict[oligo][2])
			oligostrand = oligodict[oligo][3]
			oligoage = oligodict[oligo][4]
			oligolocation = oligodict[oligo][5]
			#Skip any berglund oligos.  They don't have R values.  You will get warnings that they don't have an oligo match.
			if oligo.split('|')[0] == 'berglund':
				continue
			if chrm == oligochrm and oligostart <= kmerstart and oligostop >= kmerstop:
				match = True
				#tempoligoname = oligo.split('|')
				#roligoname = ('|').join([tempoligoname[0], tempoligoname[1], tempoligoname[2], tempoligoname[4], tempoligoname[5], tempoligoname[6], tempoligoname[7], tempoligoname[8]])
				#rvalue = rvaluedict[roligoname]
				rvalue = rvaluedict[oligo]
				if oligostrand == '+':
					outfh.write(('\t').join([oligo, oligochrm, str(oligostart), str(oligostop), oligostrand, oligoage, oligolocation, rvalue, str(kmerstart), str(kmerstop), str(branchlength), str(kmer.transcribe())]) + '\n')
				elif oligostrand == '-':
					outfh.write(('\t').join([oligo, oligochrm, str(oligostart), str(oligostop), oligostrand, oligoage, oligolocation, rvalue, str(kmerstart), str(kmerstop), str(branchlength), str(kmer.reverse_complement().transcribe())]) + '\n')
		if not match:
			print 'WARNING: no match found for {0}: {1}-{2}!'.format(chrm, kmerstart -1, kmerstop)
	
	blfh.close()
	outfh.close()

rvaluedict = getRvalues(sys.argv[4])
oligodict = getOligocoords(sys.argv[1])
addfieldstobranchlength(sys.argv[2], oligodict, rvaluedict, sys.argv[3])


