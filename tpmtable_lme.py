#Given a tpm table containing WT and KO replicate LR values, use a linear mixed effects model
#to see which genes show significant differences in WT and KO LR values.

from pandas import DataFrame as pDataFrame
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from collections import OrderedDict
import math
import argparse

pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

#define R packages
nlme = importr('nlme')
base = importr('base')
stats = importr('stats')
qv = importr('qvalue')

#define formulae
fmla = Formula('value ~ 1 + WT')
rndm = Formula('~ 1 | samples')
nullfmla = Formula('value ~ 1')
nullrndm = Formula('~1 | samples')

#Given a split line from a tpm table, convert it to a dictionary
#The indices of the line may have to be changed depending on the tpm table
def data2dict(data):
	d = {}
	d['Gene'] = [data[1]] * 6 #Gene name
	#Get tpms [WTSomaA, WTSomaB, WTSomaC, WTNeuriteA, WTNeuriteB, WTNeuriteC, KOSomaA, KOSomaB, KOSomaC, KONeuriteA, KONeuriteB, KONeuriteC]
	tpms = [float(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[8]), float(data[9]), float(data[10]), float(data[11]), float(data[12]), float(data[13])]
	#See which tpms fall below threshhold
	#If any tpm does, return none
	failfilters = [tpm < 5 for tpm in tpms] #[False, False, True, False, etc.]
	if sum(failfilters) >= 1:
		return None
	
	d['variable'] = ['WTALR', 'WTBLR', 'WTCLR', 'KOALR', 'KOBLR', 'KOCLR']
	d['value'] = [float(data[18]), float(data[19]), float(data[20]), float(data[21]), float(data[22]), float(data[23])] #LR values
	d['WT'] = [1, 1, 1, 0, 0, 0]
	d['KO'] = [0, 0, 0, 1, 1, 1]
	d['samples'] = [1, 2, 3, 4, 5, 6]
	return d

#Given a tpm line, return the pvalue from the LME
def getlmep(data):
	datadict = data2dict(data)
	#if datadict is none, that means some TPM failed length filter
	if datadict == None:
		return None
	df = pDataFrame.from_dict(datadict)
	lm_alt = nlme.lme(fmla, random = rndm, data = df, method = 'ML') #test
	lm_null = nlme.lme(nullfmla, random = nullrndm, data = df, method = 'ML') #control
	logratio = (stats.logLik(lm_alt)[0] - stats.logLik(lm_null)[0]) * 2
	pvalue = stats.pchisq(logratio, df = 1, lower_tail = False)[0]
	#format decimal
	pvalue = float('{:.2e}'.format(pvalue))
	return pvalue

#Given a list of pvalues, perform multiple hypothesis correction and return qvalues
def getqvalues(pvalues):
	#Turn list of pvalues into R vector
	pvec = FloatVector(pvalues)
	#Get qvalues object
	qobj = qv.qvalue(p = pvec)
	#qvalues are index 2 of qvalue object
	qvalues = list(qobj[2])
	#format decimal
	qvalues = [float('{:.2e}'.format(qvalue)) for qvalue in qvalues]
	return qvalues

def iteratetpmtable(tpmtable):
	counter = 0
	numberofgenes = 0
	#Get number of genes in tpmtable
	with open(tpmtable, 'r') as infh:
		for line in infh:
			numberofgenes +=1
	numberofgenes = numberofgenes - 1 #header line doesn't count

	#Using ordereddict so that matching qvalue to pvalue is easily done based on position of key
	pvaluedict = OrderedDict() #{genename : pvalue}
	with open(tpmtable, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			#Skip header
			if line[0] == 'ensembl_gene_id':
				continue
			counter +=1
			if counter % 1000 == 0:
				print 'Gene {0} of {1}...'.format(counter, numberofgenes)
			gene = line[1]
			pvalue = getlmep(line)
			#If pvalue is None, one of the tpms must have
			if pvalue == None:
				continue
			pvaluedict[gene] = pvalue

	genes = pvaluedict.keys()
	pvalues = pvaluedict.values()
	qvalues = getqvalues(pvalues)
	qvaluedict = dict(zip(genes, qvalues))
	return pvaluedict, qvaluedict

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--tpmtable', type = str, help = 'TPM table.')
	parser.add_argument('--outfile', type = str, help = 'TPM table with p and q values added.')
	args = parser.parse_args()

	pvaluedict, qvaluedict = iteratetpmtable(args.tpmtable)

	with open(args.tpmtable, 'r') as infh, open(args.outfile, 'w') as outfh:
		for line in infh:
			line = line.strip().split('\t')
			#Write header line
			if line[0] == 'ensembl_gene_id':
				line += ['deltaLR_p', 'deltaLR_q']
				outfh.write(('\t').join(line) + '\n')
				continue
			gene = line[1]
			#Some genes may have been filtered out because they did not meet tpm threshold for one or more replicates
			if gene not in pvaluedict:
				continue
			pvalue = str(pvaluedict[gene])
			qvalue = str(qvaluedict[gene])
			line += [pvalue, qvalue]
			outfh.write(('\t').join(line) + '\n')









