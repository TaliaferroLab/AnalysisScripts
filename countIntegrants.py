#Given a fastq file of integrated N15 randomers, go through, slice out the random part, and count how many unique ones we have.
#If an N15 has a hamming distance of 0 or 1 to any other already observed N15, it doesn't count as a new N15.

from collections import defaultdict
import sys
import gzip
import random

def hamming_distance(seq1, seq2):
	if len(seq1) != len(seq2):
		raise ValueError('Sequences must be same length to calculate Hamming distance')
	return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


def countrandomers(fastq):
	upstream = 'TCGCCGTGTAAGTTT'
	downstream = 'AAACCATTCCAAGTA'
	randomers = defaultdict(int) #{randomer : count}
	randomerlist = []
	i = 0
	linecounter = 0
	validreadcounter = 0
	with gzip.open(fastq) as infh:
		for line in infh:
			linecounter +=1
			if linecounter % 4000000 == 0:
				print('Analyzing read {0}...'.format(linecounter / 4))
			line = line.strip()
			i +=1
			if i == 2:
				seq = str(line)
				try:
					randomer = seq.split(upstream)[1].split(downstream)[0]
				except IndexError:
					continue
				if len(randomer) != 15:
					continue
				elif len(randomer) == 15:
					randomers[randomer] +=1
					randomerlist.append(randomer)
					validreadcounter +=1

			if i == 4:
				i = 0

	#print(randomers)
	print(len(randomers))
	print(validreadcounter)

	#return randomers, randomers.keys()
	return randomerlist, randomers

def subsamplerandomers(randomerlist, outfile):
	with open(outfile, 'w') as outfh:
		outfh.write(('\t').join(['frac', 'sampledreads', 'uniquereads']) + '\n')
		for i in range(1, 11):
			frac = i / 10
			numreads = round(len(randomerlist) * frac)
			chosenreads = random.sample(randomerlist, numreads)
			uniquereads = list(set(chosenreads))
			outfh.write(('\t').join([str(frac), str(numreads), str(len(uniquereads))]) + '\n')



def collapseseqs(randomers):
	collapsedseqs = {} #{randomer : count}
	allrandomers = randomers.keys() # a list of all of the keys in randomers
	similarseqs = [] #list of sequences that are similar to any other seq
	print('Identified {0} unique randomers.').format(len(randomers))
	x = 0
	for randomer in randomers:
		x +=1
		if x % 100 == 0:
			print('Collapsing randomer {0} of {1}.'.format(x, len(randomers)))
		hds = []
		for otherrandomer in allrandomers:
			hd = hamming_distance(randomer, otherrandomer)
			hds.append(hd)
		similarseqindicies = indices(hds, 1)
		similarcounts = 0
		for i in similarseqindicies:
			print(randomer, similarseqindicies, allrandomers[i])
			similarrandomer = allrandomers[i]
			similarseqs.append(similarrandomer)
			j = randomers[similarrandomer]
			similarcounts +=j

		if not similarseqindicies:
			collapsedseqs[randomer] = randomers[randomer]
		else:
			collapsedseqs[randomer] = randomers[randomer] + j

	similarseqs = list(set(similarseqs))

	for seq in similarseqs:
		collapsedseqs.pop(seq, None)

	print(collapsedseqs)

	return collapsedseqs

			
randomerlist = countrandomers(sys.argv[1])
subsamplerandomers(randomerlist, sys.argv[2])
#randomers_cad = countrandomers(sys.argv[1])[1]
#randomers_n2a = countrandomers(sys.argv[2])[1]
#u = set.intersection(set(randomers_cad), set(randomers_n2a))
#print(len(randomers_cad), len(randomers_n2a), len(u))
#collapseseqs(randomers)