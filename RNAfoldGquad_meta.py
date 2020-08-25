import math
import os
import argparse
#Given a gquadout file, like AlldeltaLRGenes.longest3UTR.gquadout.txt, turn the locations of the quadruplexed g's into a meta output

def parsegquadout(gquadout, numberofbins, outfile, region, txclass):
	metadict = {} #{bin : [number of times this bin existed, number of times this bin had a gquad in it]}
	#Populate metadict
	for i in range(numberofbins):
		metadict[i + 1] = [0, 0]

	with open(gquadout, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if line[0] == 'seqname':
				continue
			seqname = line[0]
			seqlen = int(line[1])

			#See which bins have the possibility to have something in them for this sequence
			activebins = []
			for i in range(seqlen):
				rawbin = ((i + 1) / float(seqlen)) * numberofbins
				smoothbin = int(math.ceil(rawbin))
				if smoothbin not in activebins:
					activebins.append(smoothbin)

			#See which bins have gquad g's in them
			gquadbins = []
			if line[4] != 'none':
				gquadpos = line[4].split(',')
				gquadpos = [int(pos) for pos in gquadpos]
				#These are 0-based. Turn them into 1-based.
				gquadpos = [pos+1 for pos in gquadpos]

				for pos in gquadpos:
					rawbin = (pos / float(seqlen)) * numberofbins
					smoothbin = int(math.ceil(rawbin))
					if smoothbin == 0:
						print seqname, pos, seqlen
					gquadbins.append(smoothbin)

			#Now add these values to metadict
			#First the bins that are present in this sequence
			for activebin in activebins:
				metadict[activebin][0] +=1

			#Now the bins that had gquads in them
			for gquadbin in gquadbins:
				metadict[gquadbin][1] +=1
				
	#Output frequencies
	if os.path.isfile(outfile) == False:
		with open(outfile, 'a') as f:
			f.write(('\t').join(['region', 'class', 'bin', 'freq']) + '\n')

	with open(outfile, 'a') as f:
		for b in metadict:
			freq = metadict[b][1] / float(metadict[b][0])
			f.write(('\t').join([region, txclass, str(b), str(freq)]) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gquadoutfile', type = str, help = 'RNAfold gquad out file.')
	parser.add_argument('--numberofbins', type = int, help = 'Number of bins for metagene.')
	parser.add_argument('--outfile', type = str, help = 'Output file.')
	parser.add_argument('--region', type = str, help = 'Tx region, usually UTR5, CDS, or UTR3.')
	parser.add_argument('--txclass', type = str, help = 'Mislocalization class.  Usually AllGenes, SD1, etc.')
	args = parser.parse_args()

	parsegquadout(args.gquadoutfile, args.numberofbins, args.outfile, args.region, args.txclass)





