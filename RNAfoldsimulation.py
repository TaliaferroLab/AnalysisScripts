import sys
import subprocess
from Bio import SeqIO
import random
import argparse

def parseRNAfold(structure):
	#This function isn't actually used.  It's been replaced by parseRNAfold2.
	RNAlength = len(structure)
	leftfacing = [] #indexes of left facing parentheses
	rightfacing = [] #indexes of right facing parentheses
	pairing = {} # {upstream base : downstream base it's paired to}
	#Get positions of parentheses
	for i, x in enumerate(structure):
		if x == '(':
			rightfacing.append(i + 1) #There is no base 0
		elif x == ')':
			leftfacing.append(i + 1)

	#Pair the left and right facing parentheses
	currentrights = []
	currentlefts = []
	currentdirection = 'right'
	for i in range(RNAlength):
		if i + 1 in rightfacing and currentdirection == 'right':
			currentrights.append(i + 1)
		if i + 1 in leftfacing:
			currentdirection = 'left'
			currentlefts.append(i + 1)
		if (i + 1 in rightfacing and currentdirection == 'left') or i + 1 == RNAlength:
			#Pair the current parenthesis
			if len(currentrights) != len(currentlefts):
				print 'ERROR!!! Number of left base pairs does not equal number of right base pairs!'
			print currentrights, currentlefts
			for j in range(len(currentrights)):
				leftbase = currentrights[j]
				rightbase = currentlefts[-1 - j]
				pairing[leftbase] = rightbase
				pairing[rightbase] = leftbase

			currentdirection = 'right' #Reset direction
			currentrights = [i + 1]
			currentlefts = []
			i = i - 1 #We will need to redo this in in the next round of the for loop

	#Any base that is not paired is now paired to 0
	for i in range(RNAlength):
		if i + 1 not in pairing:
			pairing[i + 1] = 0

	print pairing
	return pairing

def parseRNAfold2(structure):
	#Switchpoints are when you switch from being a leftstem ('(') to a right stem (')')
	#This function can handle a structure with multiple stems, which parseRNAfold cannot.
	switchpoints = [0] #always have a switchpoint at the very beginning
	stemleft = []
	stemright = []
	pairing = {}

	for i, x in enumerate(structure):
		if x == '(':
			stemleft.append(i + 1) #There is no base 0
		elif x == ')':
			stemright.append(i + 1)

	stemside = 'left'
	for i, x in enumerate(structure):
		if x == ')' and stemside == 'left':
			switchpoints.append(i) #This is the (1-based) base just upstream of the first ')' (i starts at 0)
			stemside = 'right'
		if x == '(':
			stemside = 'left'

	switchpoints.append(len(structure) + 1) #always have a switchpoint be after the last base

	for ind, switchpoint in enumerate(switchpoints):
		if switchpoint == 0:
			continue
		if switchpoint > len(structure):
			break
		leftstemsbeforeswitchpoint = []
		rightstemsafterswitchpoint = []

		#Get right stems after switchpoint but before next switchpoint
		for i in range(switchpoint, switchpoints[ind + 1]):
			if i in stemright:
				rightstemsafterswitchpoint.append(i)

		#Get left stems before switchpoint, including any leftovers
		for j in range(1, switchpoint):
			if j in stemleft:
				leftstemsbeforeswitchpoint.append(j)

		if len(leftstemsbeforeswitchpoint) > len(rightstemsafterswitchpoint):
			for k in range(len(rightstemsafterswitchpoint)):
				rightbase = rightstemsafterswitchpoint[k]
				leftbase = leftstemsbeforeswitchpoint[-1 - k]
				pairing[rightbase] = leftbase
				pairing[leftbase] = rightbase
				stemleft.remove(leftbase)
				stemright.remove(rightbase)

		elif len(leftstemsbeforeswitchpoint) < len(rightstemsafterswitchpoint):
			for n in range(len(leftstemsbeforeswitchpoint)):
				rightbase = rightstemsafterswitchpoint[n]
				leftbase = leftstemsbeforeswitchpoint[-1 - n]
				pairing[rightbase] = leftbase
				pairing[leftbase] = rightbase
				stemleft.remove(leftbase)
				stemright.remove(rightbase)

		elif len(leftstemsbeforeswitchpoint) == len(rightstemsafterswitchpoint):
			for p in range(len(leftstemsbeforeswitchpoint)):
				rightbase = rightstemsafterswitchpoint[p]
				leftbase = leftstemsbeforeswitchpoint[-1 - p]
				pairing[rightbase] = leftbase
				pairing[leftbase] = rightbase
				stemleft.remove(leftbase)
				stemright.remove(rightbase)

	#Any base that is not paired is now paired to 0
	for i in range(len(structure)):
		if i + 1 not in pairing:
			pairing[i + 1] = 0

	return pairing

			

def comparepairing(structureA, structureB):
	#Compare to "pairing" dicts to see how many base pairs are different
	#Different means they are paired to some other base, or are paired/unpaired if they were unpaired/paired before
	diff = {} # {base : 1 if different, 0 if not}
	for base in structureA:
		partnerA = structureA[base]
		partnerB = structureB[base]
		if partnerA == partnerB:
			diff[base] = 0
		elif partnerA != partnerB:
			diff[base] = 1

	return diff

def comparestructuretoprobs(structureA, probsB):
	#Given a pairing dict from a MFE structure (i.e. the MFE of the original sequence)
	diff = {} #{base : 1 if different, 0 if not}
	for base in structureA:
		partnerA = structureA[base]
		partnerB = probsB[base]
		if partnerA == 0:
			if not partnerB:
				diff[base] = 0
			elif partnerB:
				diff[base] = 1
		elif partnerA != 0:
			if partnerA in partnerB:
				diff[base] = 0
			elif partnerA not in partnerB:
				diff[base] = 1

	return diff

def updatediff(currentdiff, difftoadd):
	#Update the diff dictionary with new values
	updateddiff = {}
	for base in currentdiff:
		currentvalue = currentdiff[base]
		valuetoadd = difftoadd[base]
		updateddiff[base] = currentvalue + valuetoadd

	return updateddiff


def getstructure(seq):
	#Given a sequence, return it's structure in parentheses format
	command = 'RNAfold'
	job = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	job.stdin.write(seq)
	output = job.communicate()
	return str(output[0].split('\n')[1].split(' ')[0])

def getbpprobs(seq):
	#Given a sequence, make bp probabilities and return every pair that has a probability of at least <filter>
	pairs = {} # {base : [list of bp it is paired to with at least a minimum probability]}
	probfilter = 0.5
	#Populate dictionary
	for i in range(len(seq)):
		pairs[i + 1] = []
	
	command = 'RNAfold -p'
	job = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	job.stdin.write(seq)
	output = job.communicate()
	structure = str(output[0].split('\n')[1].split(' ')[0])
	print structure
	bpfh = open('dot.ps', 'r')
	#lines containing bp probs have 'ubox' in the 4th field
	for line in bpfh:
		line = line.strip().split(' ')
		if len(line) != 4:
			continue
		if line[3] == 'ubox':
			leftbase = int(line[0])
			rightbase = int(line[1])
			bpprob = float(line[2])
			if bpprob >= probfilter:
				pairs[leftbase].append(rightbase)
				pairs[rightbase].append(leftbase)

	bpfh.close()
	return pairs



def makeoriginalstructure(originalfasta):
	#Make the diff dict of the sequence you want to compare everything to
	for record in SeqIO.parse(originalfasta, 'fasta'):
		seq = str(record.seq)
	structure = getstructure(seq)
	originalpairs = parseRNAfold2(structure)
	return originalpairs

def compareseqs(fasta, seqlengths, originalpairs, outputfile):
	#Fasta = fasta of all mutagenized sequences
	#seqlengths = length of all fasta sequences (must all be the same)
	#originalpairs = structure of your standard seq to compare to (from makeoriginalstructure)
	#First make the starting diff dictionary
	startingdiff = {}
	for i in range(seqlengths):
		startingdiff[i + 1] = 0

	counter = 0

	#Now go through each sequence in the fasta 
	for record in SeqIO.parse(fasta, 'fasta'):
		print 'Folding {0}...'.format(record.id)
		counter += 1
		if counter % 1000 == 0:
			print 'Folding sequence {0}...'.format(counter) 
		seq = str(record.seq.transcribe())
		structure = getstructure(seq)
		print record.id, structure
		pairs = parseRNAfold2(structure)
		difftoadd = comparepairing(originalpairs, pairs)
		if counter == 1:
			updateddiff = updatediff(startingdiff, difftoadd)
		elif counter >1:
			updateddiff = updatediff(updateddiff, difftoadd)

	outfh = open(outputfile, 'w')
	outfh.write('base' + '\t' + 'structure_changes' + '\t' + 'original' + '\n')
	for base in sorted(updateddiff):
		if originalpairs[base] == 0: #unpaired in original
			original = '.'
		elif originalpairs[base] != 0: #paired in original
			original = '|'
		outfh.write(str(base) + '\t' + str(updateddiff[base]) + '\t' + original + '\n')
	outfh.close()

def compareseqs_probs(fasta, seqlengths, originalpairs, outputfile):
	#Fasta = fasta of all mutagenized sequences
	#seqlengths = length of all fasta sequences (must all be the same)
	#originalpairs = structure of your standard seq to compare to (from makeoriginalstructure)
	#First make the starting diff dictionary
	startingdiff = {}
	for i in range(seqlengths):
		startingdiff[i + 1] = 0

	counter = 0

	#Not go through each sequence in the fasta
	for record in SeqIO.parse(fasta, 'fasta'):
		print 'Folding {0}...'.format(record.id)
		counter += 1
		if counter % 1000 == 0:
			print 'Folding sequence {0}...'.format(counter)
		seq = str(record.seq.transcribe())
		probs = getbpprobs(seq)
		difftoadd = comparestructuretoprobs(originalpairs, probs)
		if counter == 1:
			updateddiff = updatediff(startingdiff, difftoadd)
		elif counter > 1:
			updateddiff = updatediff(updateddiff, difftoadd)

	outfh = open(outputfile, 'w')
	outfh.write('base' + '\t' + 'structure_changes' + '\t' + 'original' + '\n')
	for base in sorted(updateddiff):
		if originalpairs[base] == 0: #unpaired in original
			original = '.'
		elif originalpairs[base] != 0: #paired in original
			original = '|'
		outfh.write(str(base) + '\t' + str(updateddiff[base]) + '\t' + original + '\n')
	outfh.close()


def makemutatedfasta(inputfasta, iterations):
	mutdict = {'A' : ['C','U','G'], 'G' : ['A','U','C'], 'C' : ['A','U','G'], 'U' : ['A','C','G'], 'T' : ['A','C','G']}
	lottery = []
	outfh = open('MutatedSequences.fasta', 'w')
	for i in range(1000):
		lottery.append(i + 1)

	#Get the original sequence
	for record in SeqIO.parse(inputfasta, 'fasta'):
		seq = str(record.seq)

	for i in range(int(iterations)):
		if (i + 1) % 10000 == 0:
			print 'Creating mutated sequence {0} of {1}...'.format(i + 1, iterations)
		seqtitle = 'Seq' + str(i + 1)
		currentseq = ''
		for nt in seq:
			ticket = random.choice(lottery)
			if ticket <= 949:
				currentseq += nt
			elif ticket > 949 and ticket <= 966:
				mutation = mutdict[nt][0]
				currentseq += mutation
			elif ticket > 966 and ticket <= 983:
				mutation = mutdict[nt][1]
				currentseq += mutation
			elif ticket > 983 and ticket <= 1000:
				mutation = mutdict[nt][2]
				currentseq += mutation

		outfh.write('>' + seqtitle + '\n' + currentseq + '\n')

	print 'Done making mutated fasta!!'

	outfh.close()



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--originalsequence', type = str, help = 'Original sequence in fasta format.')
	parser.add_argument('--iterations', type = int, help = 'Number of mutated sequences to make.')
	parser.add_argument('--mode', choices = ['MFE','prob'], type = str, help = 'Compare MFE to MFE or MFE of original sequence to probs in mutated sequences?')
	parser.add_argument('--output', type = str, help = 'Output file of structure changes.')
	args = parser.parse_args()

	#Get length of original sequence
	for record in SeqIO.parse(args.originalsequence, 'fasta'):
		seqlength = len(str(record.seq))

	#Make mutated fasta
	print 'Making mutated sequences...'
	makemutatedfasta(args.originalsequence, args.iterations)

	#Define basepairs in original fasta
	print 'Defining basepairs in original sequence...'
	originalpairs = makeoriginalstructure(args.originalsequence)

	if args.mode == 'MFE':
		#Fold and compare mutated sequences
		compareseqs('MutatedSequences.fasta', int(seqlength), originalpairs, args.output)

	elif args.mode == 'prob':
		compareseqs_probs('MutatedSequences.fasta', int(seqlength), originalpairs, args.output)

