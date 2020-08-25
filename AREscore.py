#http://arescore.dkfz.de/info.html

from Bio import SeqIO
import re

motif = 'ATTTA'

def getmatches(motif, seq):
	matches = []
	for match in re.finditer(motif, seq):
		matches.append([match.start(), match.end()])

def getscore(matches):
	score = 0
	for match in matches:
		score += 1
	if len(matches) > 1:
		for ind, match in enumerate(matches):
			if ind <= len(matches) - 1:
				leftmatchend = matches[ind][1]
				rightmatchstart = matches[ind + 1][0]
				intermatchdistance = rightmatchstart - leftmatchend


with open(fasta, 'r') as infh:
	for record in SeqIO.parse(fasta, 'fasta'):
		seq = str(record.seq)
		matches = getmatches(motif, seq)

