#Given a modelgff, subsample another gff so that the distribution of lengths matches the distribution of modelgff
import gffutils
import random
import math
import os
import sys

def subsampleLength(modelgff, subsamplegff, numberofbins, subsamplecheckoutput):
	modelbins = {}
	subsamplebins = {}
	modelgfffeatures = 0
	subsamplegfffeatures = 0
	subsampledfeatures = []
	#Make gff databases

	#Gff you want to model after
	print 'Indexing gff...'
	gff_fn = modelgff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	modeldb = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Gff to be subsampled
	print 'Indexing gff...'
	gff_fn = subsamplegff
	db_fn = os.path.basename(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, verbose = True)
	subsampledb = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	#Populate modelbins

	#First get min and max length of modelgff
	featurelengths = []
	genes = modeldb.features_of_type('gene')
	for gene in genes:
		featurelength = 0
		for exon in modeldb.children(gene, featuretype = 'exon'):
			exonlength = exon.end - exon.start
			featurelength += exonlength

		featurelengths.append(featurelength)

	minlength, maxlength = min(featurelengths), max(featurelengths)

	#Now populate the bins
	genes = modeldb.features_of_type('gene')
	for gene in genes:
		modelgfffeatures +=1
		featurelength = 0
		for exon in modeldb.children(gene, featuretype = 'exon'):
			exonlength = exon.end - exon.start
			featurelength += exonlength

		lowerbound = minlength
		upperbound = lowerbound + ((maxlength - minlength) / numberofbins)
		while upperbound <= maxlength:
			if featurelength >= lowerbound and featurelength < upperbound:
				if '{0}_to_{1}'.format(lowerbound, upperbound) in modelbins:
					modelbins['{0}_to_{1}'.format(lowerbound, upperbound)].append(str(gene.id))
					break
				else:
					modelbins['{0}_to_{1}'.format(lowerbound, upperbound)] = [str(gene.id)]
					break

			else:
				#Move bounds up by one bin
				lowerbound += (maxlength - minlength) / numberofbins
				upperbound += (maxlength - minlength) / numberofbins

	#Populate subsamplebins
	genes = subsampledb.features_of_type('gene')
	for gene in genes:
		featurelength = 0
		subsamplegfffeatures +=1
		for exon in subsampledb.children(gene, featuretype = 'exon'):
			exonlength = exon.end - exon.start
			featurelength += exonlength
		
		lowerbound = minlength
		upperbound = lowerbound + ((maxlength - minlength) / numberofbins)
		while upperbound <= maxlength:
			if featurelength >= lowerbound and featurelength < upperbound:
				if '{0}_to_{1}'.format(lowerbound, upperbound) in subsamplebins:
					subsamplebins['{0}_to_{1}'.format(lowerbound, upperbound)].append(str(gene.id))
					break
				else:
					subsamplebins['{0}_to_{1}'.format(lowerbound, upperbound)] = [str(gene.id)]
					break

			else:
				#Move bound up by one bin
				lowerbound += (maxlength - minlength) / numberofbins
				upperbound += (maxlength - minlength) / numberofbins

	#Number of features in each gff...used for calculating density
	modelgfffeatures = float(modelgfffeatures)
	subsamplegfffeatures = float(subsamplegfffeatures)

	ratios = []
	for modelbin in modelbins:
		modelbinpop = float(len(modelbins[modelbin]))
		subsamplebinpop = float(len(subsamplebins[modelbin]))
		#What's the relative population of this bin in the subsample and model?
		ratio = int(math.floor(subsamplebinpop / modelbinpop))
		ratios.append(ratio)
	#We want the smallest relative population ratio because that's what we are going to 
	#multiply all the modelbin populations by
	binfactor = min(ratios)

	for modelbin in modelbins:
		modelbinpop = float(len(modelbins[modelbin]))
		subsample_features_to_pick = int(modelbinpop * binfactor)
		if modelbin in subsamplebins:
			random_subsampled_features = random.sample(subsamplebins[modelbin], subsample_features_to_pick)
			subsampledfeatures += random_subsampled_features


	print 'There were {0} features in the model gff and {1} in the gff to be subsampled. {2} features were chosen in the sampling.'.format(modelgfffeatures, subsamplegfffeatures, len(subsampledfeatures))

	'''
	#Write subsamplecheckoutput so that you can verify that the length distributions are the same
	with open(subsamplecheckoutput, 'w') as f:
		f.write('sample' + '\t' + 'ID' + '\t' + 'length' + '\n')
		genes = modeldb.features_of_type('gene')
		for gene in genes:
			ID = str(gene.id)
			featurelength = 0
			for exon in modeldb.children(gene, featuretype = 'exon'):
				featurelength += exon.end - exon.start

			f.write('model' + '\t' + ID + '\t' + str(featurelength) + '\n')
		
		genes = subsampledb.features_of_type('gene')
		for gene in genes:
			ID = str(gene.id)
			featurelength = 0
			for exon in subsampledb.children(gene, featuretype = 'exon'):
				featurelength += exon.end - exon.start
			f.write('old' + '\t' + ID + '\t' + str(featurelength) + '\n')
			if ID in subsampledfeatures:
				f.write('subsampled' + '\t' + ID + '\t' + str(featurelength) + '\n')
	'''

	#Return list of features that from subsamplegff that were chosen
	return subsampledfeatures



if __name__ == '__main__':
	subsampleLength(sys.argv[1], sys.argv[2], 40, 'test.txt')





