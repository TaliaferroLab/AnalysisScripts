#python3

#This scripts takes a directory full of salmon output directories and outputs a file of tpms where rows are transcripts and columns
#are samples for every sample in the directory.

import os
import pandas as pd
import sys
from functools import reduce

def join_dfs(ldf, rdf):
	return ldf.join(rdf, how = 'inner')

def salmon2suppa(salmonoutdir):
	dfs = []
	for d in os.listdir(salmonoutdir):
		tpmfile = os.path.join(os.path.abspath(salmonoutdir), d, 'quant.sf')
		tpmdf = pd.read_csv(tpmfile, sep = '\t', header = 0, index_col = 'Name')
		del tpmdf.index.name
		tpmdf = tpmdf.loc[:,['TPM']] #this has to be a list so you get a dataframe, not a series
		tpmdf = tpmdf.rename(columns = {'TPM' : d})

		dfs.append(tpmdf)

	#Merge dataframes
	alltpmdfs = reduce(join_dfs, dfs)
	
	outpath = os.path.join(os.path.abspath(salmonoutdir), 'tpmsforsuppa.txt')
	alltpmdfs.to_csv(outpath, sep = '\t', float_format = '%g')

	#Remove tab space in header over transcript ids
	with open(outpath, 'r') as infh, open('tabremoved.txt', 'w') as outfh:
		for line in infh:
			line = line.strip().split('\t')
			outfh.write(('\t').join(line) + '\n')

	os.rename(os.path.join(os.path.abspath(salmonoutdir), 'tabremoved.txt'), 'tpmsforsuppa.txt')

salmon2suppa(sys.argv[1])