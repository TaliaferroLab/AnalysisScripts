#python3
import pandas as pd
import seaborn as sns

def getbpprobs(bpprobfile):
	#probabilities of all possible interactions
	bpprobs = {} #nested dictionary {position1 : {position2 : prob, position3: prob}, position2: {position1 : prob, position3 : prob}}
	with open(bpprobfile, 'r') as infh:
		for line in infh:
			line = line.strip().split(' ')
			if len(line) != 4:
				continue
			if line[3] == 'ubox':
				leftbase = int(line[0])
				rightbase = int(line[1])
				bpprob = float(line[2])**2
				if leftbase not in bpprobs:
					bpprobs[leftbase] = {}
				if rightbase not in bpprobs:
					bpprobs[rightbase] = {}
				bpprobs[leftbase][rightbase] = bpprob
				bpprobs[rightbase][leftbase] = bpprob

	#For each position, if a given position is not listed in its dict, add it and set probability to 0
	positions = list(bpprobs.keys())
	for positionx in positions:
		for positiony in positions:
			if positiony not in bpprobs[positionx]:
				bpprobs[positionx][positiony] = 0.0

	#Turn nested dictionary into a list
	bpprobs_lists = {} #{pos : [list of probs, IN ORDER]}
	for posx in bpprobs:
		bpprobs_lists[posx] = []
		for posy in sorted(bpprobs[posx].keys()):
			bpprobs_lists[posx].append(bpprobs[posx][posy])

	#Get summed bp prob for a base across all other bases
	bpprobs_sum = {} #{position : sum of all bp probs}
	for pos in bpprobs_lists.keys():
		probs = bpprobs_lists[pos]
		sumprob = sum(probs)
		bpprobs_sum[pos] = sumprob

	return bpprobs_lists, bpprobs_sum



getbpprobs('7skbpprob.txt')