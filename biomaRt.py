#Use rpy2 and biomaRt to get transcript/gene relationships

from pandas import DataFrame as pDataFrame
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

#import r package
bm = importr('biomaRt')


#Get relationships between ensembl gene IDs and mouse common gene names
def getens2common():
	#get mart
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	t2g = bm.getBM(attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'refseq_mrna'], mart = mart)

	#convert t2g from an R dataframe to a pandas dataframe
	t2g = pandas2ri.ri2py(t2g)

	#zip ensembl gene name and common gene name together
	ens2common = dict(zip(t2g['ensembl_gene_id'], t2g['external_gene_name']))

	i = 0
	for gene in ens2common:
		if ens2common[gene] != '':
			i +=1

	print 'Found names for {0} of {1} ensembl gene identifiers.'.format(i, len(ens2common))

	return ens2common

#Get relationship between mouse ensembl geneIDs and human common gene names
def getens2humancommon():
	#get mart
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	t2g = bm.getBM(attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'hsapiens_homolog_associated_gene_name'], mart = mart)

	#convert t2g from an R dataframe to a pandas dataframe
	t2g = pandas2ri.ri2py(t2g)

	#zip ensembl gene name and common gene name together
	ens2humancommon = dict(zip(t2g['ensembl_gene_id'], t2g['hsapiens_homolog_associated_gene_name']))

	i = 0
	for gene in ens2humancommon:
		if ens2humancommon[gene] != '':
			i +=1

	print 'Found human gene names for {0} of {1} ensembl gene identifiers.'.format(i, len(ens2humancommon))

	return ens2humancommon


#Get relationship between human ensembl geneIDs and human common gene names
def gethens2humancommon():
	#get mart
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = 'www.ensembl.org')
	t2g = bm.getBM(attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'], mart = mart)

	#convert t2g from an R dataframe to a pandas dataframe
	t2g = pandas2ri.ri2py(t2g)

	#zip ensembl gene name and common gene name together
	hens2humancommon = dict(zip(t2g['ensembl_gene_id'], t2g['external_gene_name']))

	i = 0
	for gene in hens2humancommon:
		if hens2humancommon[gene] != '':
			i +=1

	print 'Found human gene names for {0} of {1} ensembl gene identifiers.'.format(i, len(hens2humancommon))

	return hens2humancommon

#Get relationships between human refseq transcript IDs and ensembl gene IDs
def gethrefseq2humanens():
	#get mart
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = 'www.ensembl.org')
	t2g = bm.getBM(attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'refseq_mrna'], mart = mart)

	#convert t2g from an R dataframe to a pandas dataframe
	t2g = pandas2ri.ri2py(t2g)

	#zip ensembl gene name and common gene name together
	hrefseq2humanens = dict(zip(t2g['refseq_mrna'], t2g['ensembl_gene_id']))

	i = 0
	for mRNA in hrefseq2humanens:
		if hrefseq2humanens[mRNA] != '':
			i +=1

	print 'Found human ensembl gene names for {0} of {1} refseq mrna identifiers.'.format(i, len(hrefseq2humanens))

	return hrefseq2humanens

def gethenst2hensgene():
	#get mart
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = 'www.ensembl.org')
	t2g = bm.getBM(attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'refseq_mrna'], mart = mart)

	#convert t2g from an R dataframe to a pandas dataframe
	t2g = pandas2ri.ri2py(t2g)

	#zip ensembl gene name and common gene name together
	henst2hensgene = dict(zip(t2g['ensembl_transcript_id'], t2g['ensembl_gene_id']))

	i = 0
	for txname in henst2hensgene:
		if henst2hensgene[txname] != '':
			i +=1

	print 'Found human ensembl gene names for {0} of {1} human ensembl mrna identifiers.'.format(i, len(henst2hensgene))

	return henst2hensgene

#Get relationships between mouse ensembl gene ids and human ensembl gene ids
def getmouseens2humanens():
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	m2h = bm.getBM(attributes = ['ensembl_gene_id', 'hsapiens_homolog_ensembl_gene'], mart = mart)

	#convert hm2h from an R dataframe to a pandas DataFrame
	m2h = pandas2ri.ri2py(m2h)

	#zip mouse ensembl ID and human homolog ensembl ID together
	mouseens2humanens = dict(zip(m2h['ensembl_gene_id'], m2h['hsapiens_homolog_ensembl_gene']))

	i = 0
	for mouseens in mouseens2humanens:
		if mouseens2humanens[mouseens] != '':
			i +=1

	print 'Found human homologs for {0} of {1} mouse ensembl gene ids.'.format(i, len(mouseens2humanens))

	return mouseens2humanens

#Get relationships between mouse ensembl gene ids and human ensembl gene ids
def getmouseens2dmelens():
	mart = bm.useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host = 'www.ensembl.org')
	m2d = bm.getBM(attributes = ['ensembl_gene_id', 'dmelanogaster_homolog_ensembl_gene'], mart = mart)

	#convert hm2h from an R dataframe to a pandas DataFrame
	m2d = pandas2ri.ri2py(m2d)

	#zip mouse ensembl ID and human homolog ensembl ID together
	mouseens2dmelens = dict(zip(m2d['ensembl_gene_id'], m2d['dmelanogaster_homolog_ensembl_gene']))

	i = 0
	for mouseens in mouseens2dmelens:
		if mouseens2dmelens[mouseens] != '':
			i +=1

	print 'Found Drosophila homologs for {0} of {1} mouse ensembl gene ids.'.format(i, len(mouseens2dmelens))

	return mouseens2dmelens



