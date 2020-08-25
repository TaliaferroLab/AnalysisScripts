#Given a PSI table (produced by MISOSigTimecourse_v2.0.py), calculate how many sig events there are for each sample
#in each event type, and what fraction of those events have positive delta psi values.  Delta PSI Filters for "sig event" 
#can be modified (usually 0 or 0.1). For each sample, only considers events that are marked as significant in the 
#'SigComparisons' column output by MISOSigTimecourse_v2.0.py.  This is usually BF >= 10 and an optional dPSI filter.

#Any sample name that has a lower case 'v' in it is going to cause problems.

#Outputs two files. One contains the number of significant events in each sample.  The other contains the fraction of those
#events that have positive delta psis.

import argparse
from numpy import median
from scipy.stats import chisquare

def parsetable(psitable):
	psidict = {} # {sample : {eventtype : {event : psivalue}}}
	samples = [] #same order as they appear in the table
	fh = open(psitable, 'r')
	for line in fh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			for sample in line[2:-1]: #Skip 'Event', 'EventType', 'SigComparisons'
				samples.append(sample)
				psidict[sample] = {}
		else:
			for i in range(len(samples)):
				sample = samples[i]
				eventtype = line[1]
				if psidict[sample].has_key(eventtype) == False:
					psidict[sample][eventtype] = {}
				event = line[0]
				psivalue = line[i + 2]
				psidict[sample][eventtype][event] = float(psivalue)

	fh.close()

	return psidict

def getfractions(psitable, psidict):
	fractiondict = {} # {sample : {eventtype : [sig events, number of events with positive delta PSI, number of events with negative delta psi]}}
	fh = open(psitable, 'r')
	for line in fh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			continue
		event = line[0]
		eventtype = line[1]
		sigcomparisons = line[-1].split(';')
		for sigcomparison in sigcomparisons: #Assuming that these are written as ExperimentalvControl
			if sigcomparison == 'None':
				continue
			sample1 = sigcomparison.split('v')[0] #Experimental sample....any sample name with a lower case v is going to be a problem
			sample2 = sigcomparison.split('v')[1] #Control sample
			psivalue1 = psidict[sample1][eventtype][event]
			psivalue2 = psidict[sample2][eventtype][event]
			if psivalue1 and psivalue2:
				if fractiondict.has_key(sample1) == False:
					fractiondict[sample1] = {}
				if fractiondict[sample1].has_key(eventtype) == False:
					fractiondict[sample1][eventtype] = [0,0,0]
			if abs(psivalue1 - psivalue2) >= 0:  #Can change this
				fractiondict[sample1][eventtype][0] +=1
			if psivalue2 - psivalue1 <= 0:  #Can change this
				fractiondict[sample1][eventtype][1] +=1
			if psivalue2 - psivalue1 >= 0:
				fractiondict[sample1][eventtype][2] +=1

	fh.close()
	return fractiondict

def getmediandeltapsis(psitable, psidict):
	#Not normally used
	mediandict = {} # {sample : {eventtype : [deltapsis of significant events]}}
	outdict = {} # {sample : {eventtype : median delta psi of significant events}}
	fh = open(psitable, 'r')
	for line in fh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			continue
		event = line[0]
		eventtype = line[1]
		sigcomparisons = line[-1].split(';')
		for sigcomparison in sigcomparisons: #Assuming that these are written as ExperimentalvControl
			sample1 = sigcomparison.split('v')[0] #Experimental sample....any sample name with a lower case v is going to be a problem
			sample2 = sigcomparison.split('v')[1] #Control sample
			psivalue1 = psidict[sample1][eventtype][event]
			psivalue2 = psidict[sample2][eventtype][event]
			if psivalue1 and psivalue2:
				if abs(psivalue1 - psivalue2) >= 0.1: #can change this
					deltapsi = float(psivalue1 - psivalue2)
					if mediandict.has_key(sample1) == False:
						mediandict[sample1] = {}
					if mediandict[sample1].has_key(eventtype) == False:
						mediandict[sample1][eventtype] = []
					mediandict[sample1][eventtype].append(deltapsi)

	fh.close()
	for sample in mediandict:
		if outdict.has_key(sample) == False:
			outdict[sample] = {}
		for eventtype in mediandict[sample]:
			if outdict[sample].has_key(eventtype) == False:
				outdict[sample][eventtype] = median(mediandict[sample][eventtype])

	outfh = open('mediantest.txt', 'w')
	outfh.write('Sample' + '\t' + 'ALE' + '\t' + 'TandemUTR' + '\n')
	for sample in outdict:
		outfh.write(('\t').join([sample, str(outdict[sample]['ALE']), str(outdict[sample]['TandemUTR'])]) + '\n')
	return outdict


def writefractions(fractiondict, outfile):
	outfh = open(outfile, 'w')
	outfh.write(('\t').join(['Sample','AFE','AFESigEvents','ALE','ALESigEvents','SE','SESigEvents','TandemUTR','TandemUTRSigEvents']) + '\n')
	eventtypes = ['AFE','ALE','SE','TandemUTR']
	for sample in sorted(fractiondict):
		if len(fractiondict[sample]) >= 4: #only consider samples with at least one sig event in each eventtype
			outfh.write(sample)
			for eventtype in eventtypes:
				fracpos = round(fractiondict[sample][eventtype][1] / float(fractiondict[sample][eventtype][1] + fractiondict[sample][eventtype][2]), 3)
				sigevents = fractiondict[sample][eventtype][0]
				outfh.write('\t' + str(fracpos) + '\t' + str(sigevents))
			outfh.write('\n')
	outfh.close()


def chisq(infile):
	#Given an output from this script, calculate chi square p values based on a 50/50 expectation of having a pos or neg delta PSI value.
	#Takes the number of sig events and fraction with pos delta psi (given in output of this script).
	#Prints Sample / EventType combinations that have chi square pvalues less than 0.05

	infh = open(infile, 'r')
	eventtypes = ['AFE','ALE','SE','TandemUTR']
	eventtypeindex = {1 : 'AFE', 3 : 'ALE', 5 : 'MXE', 7 : 'RI', 9 : 'SE', 11 : 'TandemUTR'} #these are the indexes of the fraction of events with pos delta psi for each event type
	for line in infh:
		line = line.strip().split('\t')
		if line[0] == 'Sample':
			continue
		sample = line[0]
		for i in range(1,11,2):
			eventtype = eventtypeindex[i]
			fracpos = float(line[i])
			sigevents = float(line[i+1])
			posevents = fracpos * sigevents
			negevents = (1 - fracpos) * sigevents
			pvalue = chisquare([posevents, negevents])[1]
			if pvalue > 0.05 and (eventtype == 'ALE' or eventtype == 'TandemUTR'):
				print 'The {0} events in {1} have a chi square pvalue of {2}.'.format(eventtype, sample, str(pvalue))

	infh.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--psitable', type = str, help = 'PSI table created by MISOSigTimecourse_v2.py')
    parser.add_argument('--output', type = str, help = 'Output file containing the number of significant events in each sample in each event type and the fraction that has positive delta psi values.')
    parser.add_argument('--chisquareinput', type = str, help = 'Input file that is the output of this script. Will use chi square test to determine which fractions are significantly different from 50/50.') 
    args = parser.parse_args()

    if not args.chisquareinput:
    	psidict = parsetable(args.psitable)
    	fractiondict = getfractions(args.psitable, psidict)
    	writefractions(fractiondict, args.output)

    elif args.chisquareinput:
    	chisq(args.chisquareinput)