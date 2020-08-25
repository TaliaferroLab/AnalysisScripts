import numpy
import sys


def fitline(data1, data2):
	#Fit a line for the data. Data1 is on xaxis, data2 on y axis. Return slope and intercept of that line.
	if len(data1) == 0 or len(data2) == 0:
		print 'hullo'
		return 1, 0
	coefficients = numpy.polyfit(data1, data2, 1)
	slope = coefficients[0]
	yint = coefficients[1]

	return slope, yint

def getresidual(xval, yval, slope, intercept):
	predictedyval = (xval * slope) + intercept
	residual = yval - predictedyval

	return residual

def calculateresiduals(psitable, sample1soma, sample1axon, sample2soma, sample2axon):
	#Sample1 is the 'x' value and sample2 is the 'y' value
	PSIs = {} # {eventtype : [[sample1psis], [sample2psis]]}
	fits = {} # {eventtype : [slope, yint]}
	residualcolname = 'Hr12C_v_Hr12D_residual'
	infh = open(psitable, 'r')
	header = infh.readline().strip().split('\t')
	infh.close()

	#Figure out which columns you want
	for idx, colname in enumerate(header):
		if colname == sample1soma:
			sample1somaidx = idx
		if colname == sample1axon:
			sample1axonidx = idx
		if colname == sample2soma:
			sample2somaidx = idx
		if colname == sample2axon:
			sample2axonidx = idx

	try:
		sample1somaidx, sample1axonidx, sample2somaidx, sample2axonidx
	except NameError:
		print 'Could not find sample names in PSI table!!'
		sys.exit()

	#Get the dPSI values
	infh = open(psitable, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			continue
		eventtype = line[1]
		SigComparisons = line[10]
		if eventtype not in PSIs:
			PSIs[eventtype] = [[],[]]
		sample1somapsi = float(line[sample1somaidx])
		sample1axonpsi = float(line[sample1axonidx])
		sample2somapsi = float(line[sample2somaidx])
		sample2axonpsi = float(line[sample2axonidx])
		sample1dpsi = sample1axonpsi - sample1somapsi
		sample2dpsi = sample2axonpsi - sample2somapsi
		#For the purposes of residuals, only consider events that were localized in both reps of
		#Hr0 and/or both reps of Hr12
		if ('Hr0AxonCvHr0SomaC' in SigComparisons and 'Hr0AxonDvHr0SomaD' in SigComparisons) or ('Hr12AxonCvHr12SomaC' in SigComparisons and 'Hr12AxonDvHr12SomaD' in SigComparisons):
			PSIs[eventtype][0].append(sample1dpsi)
			PSIs[eventtype][1].append(sample2dpsi)

	infh.close()

	#Fit lines for the dPSI values
	for eventtype in PSIs:
		slope, yint = fitline(PSIs[eventtype][0], PSIs[eventtype][1])
		print eventtype, slope, yint
		fits[eventtype] = [slope, yint]

	#Write a new PSI table with dPSI residuals added
	infh = open(psitable, 'r')
	outfh = open('/Users/mtaliaferro/Documents/MIT/Localization/CADDepolarization/May2015/MISO/PSITable_residuals6.txt', 'w')
	for line in infh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			outfh.write(('\t').join(line) + '\t' + residualcolname + '\n')
			continue
		eventtype = line[1]
		slope = fits[eventtype][0]
		intercept = fits[eventtype][1]
		sample1somapsi = float(line[sample1somaidx])
		sample1axonpsi = float(line[sample1axonidx])
		sample2somapsi = float(line[sample2somaidx])
		sample2axonpsi = float(line[sample2axonidx])
		sample1dpsi = sample1axonpsi - sample1somapsi
		sample2dpsi = sample2axonpsi - sample2somapsi
		residual = getresidual(sample1dpsi, sample2dpsi, slope, intercept)
		outfh.write(('\t').join(line) + '\t' + str(residual) + '\n')

	infh.close()
	outfh.close()

def getsigchanging(residualtable):
	passingevents = {} # {eventtype : [events that pass filters]}
	infh = open(residualtable, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			continue
		eventname = line[0]
		eventtype = line[1]
		SigComparisons = line[10]
		cr1 = float(line[11])
		cr2 = float(line[16])
		er1 = float(line[12])
		er2 = float(line[13])
		er3 = float(line[14])
		er4 = float(line[15])

		controlresiduals = [cr1, cr2]
		experimentalresiduals = [er1, er2, er3, er4]
		abscontrolresiduals = [abs(cr1), abs(cr2)]
		absexperimentalresiduals = [abs(er1), abs(er2), abs(er3), abs(er4)]

		if ('Hr0AxonCvHr0SomaC' in SigComparisons and 'Hr0AxonDvHr0SomaD' in SigComparisons) or ('Hr12AxonCvHr12SomaC' in SigComparisons and 'Hr12AxonDvHr12SomaD' in SigComparisons):
			siglocalized = True
		else:
			siglocalized = False

		#check if er samples all have same sign
		if all(value >= 0 for value in experimentalresiduals) or all(value <= 0 for value in experimentalresiduals) == True:
			samesign = True
		else:
			samesign = False

		if min(absexperimentalresiduals) >= max(abscontrolresiduals):
			abspass = True
		else:
			abspass = False

		if samesign == True and abspass == True and siglocalized == True:
			if passingevents.has_key(eventtype) == False:
				passingevents[eventtype] = []
			passingevents[eventtype].append(eventname)

	for eventtype in passingevents:
		print eventtype, len(passingevents[eventtype])

	outfh = open('ALEpassedevents.txt', 'w')
	for event in passingevents['ALE']:
		outfh.write(event + '\n')
	outfh.close()

	infh.close()




#calculateresiduals(sys.argv[1], 'Hr0SomaC', 'Hr0AxonC', 'Hr0SomaD', 'Hr0AxonD')
getsigchanging('PSITable_residuals.txt')



	