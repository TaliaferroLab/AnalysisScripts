import numpy as np
import sys


def filterdpsitable(ddpsitable, outfile):
	counter = 0
	passingcounter = 0
	unchangedeventcounter = 0
	outfh = open(outfile, 'w')
	infh = open(ddpsitable, 'r')
	for line in infh:
		line = line.strip().split('\t')
		if line[0] == 'Event':
			outfh.write(('\t').join(line) + '\n')
			continue

		event = line[0]
		eventtype = line[1]
		KOaxonApsi = float(line[2])
		KOaxonBpsi = float(line[3])
		KOaxonCpsi = float(line[4])
		KOsomaApsi = float(line[5])
		KOsomaBpsi = float(line[6])
		KOsomaCpsi = float(line[7])
		WTaxonApsi = float(line[8])
		WTaxonBpsi = float(line[9])
		WTaxonCpsi = float(line[10])
		WTsomaApsi = float(line[11])
		WTsomaBpsi = float(line[12])
		WTsomaCpsi = float(line[13])
		sigcomp = line[14]

		#Was this event significant (BF >=5) in these comparisons (WT v KO)?
		#Only considering soma samples
		WTintrasig = 0 #sig between wt replicates
		KOintrasig = 0 #sig between ko replicates
		intersig = 0 #sig between wt/ko replicates
		mindpsi = 0.05
		dpsisigfilter = False
		
		
		if 'WTSomaAvWTSomaB' in sigcomp and abs(WTsomaBpsi - WTsomaApsi) >= mindpsi:
			WTintrasig +=1
		if 'WTSomaAvWTSomaC' in sigcomp and abs(WTsomaCpsi - WTsomaApsi) >= mindpsi:
			WTintrasig +=1
		if 'WTSomaBvWTSomaC' in sigcomp and abs(WTsomaCpsi - WTsomaBpsi) >= mindpsi:
			WTintrasig +=1
		if 'KOSomaAvKOSomaB' in sigcomp and abs(KOsomaBpsi - KOsomaApsi) >= mindpsi:
			KOintrasig +=1
		if 'KOSomaAvKOSomaC' in sigcomp and abs(KOsomaCpsi - KOsomaApsi) >= mindpsi:
			KOintrasig +=1
		if 'KOSomaBvKOSomaC' in sigcomp and abs(KOsomaCpsi - KOsomaBpsi) >= mindpsi:
			KOintrasig +=1
		if 'WTSomaAvKOSomaA' in sigcomp and abs(KOsomaApsi - WTsomaApsi) >= mindpsi:
			intersig +=1
		if 'WTSomaAvKOSomaB' in sigcomp and abs(KOsomaBpsi - WTsomaApsi) >= mindpsi:
			intersig +=1
		if 'WTSomaAvKOSomaC' in sigcomp and abs(KOsomaCpsi - WTsomaApsi) >= mindpsi:
			intersig +=1
		if 'WTSomaBvKOSomaA' in sigcomp and abs(KOsomaApsi - WTsomaBpsi) >= mindpsi:
			intersig +=1
		if 'WTSomaBvKOSomaB' in sigcomp and abs(KOsomaBpsi - WTsomaBpsi) >= mindpsi:
			intersig +=1
		if 'WTSomaBvKOSomaC' in sigcomp and abs(KOsomaCpsi - WTsomaBpsi) >= mindpsi:
			intersig +=1
		if 'WTSomaCvKOSomaA' in sigcomp and abs(KOsomaApsi - WTsomaCpsi) >= mindpsi:
			intersig +=1
		if 'WTSomaCvKOSomaB' in sigcomp and abs(KOsomaBpsi - WTsomaCpsi) >= mindpsi:
			intersig +=1
		if 'WTSomaCvKOSomaC' in sigcomp and abs(KOsomaCpsi - WTsomaCpsi) >= mindpsi:
			intersig +=1

		if WTintrasig == 0 and KOintrasig == 0 and intersig >= 1:
			dpsisigfilter = True

		#Is the lowest PSI of one condition higher than the highest PSI of the other?
		dpsifilter = False
		WTpsis = [WTsomaApsi, WTsomaBpsi, WTsomaCpsi]
		KOpsis = [KOsomaApsi, KOsomaBpsi, KOsomaCpsi]
		if min(WTpsis) >= max(KOpsis) or max(WTpsis) <= min(KOpsis):
			dpsifilter = True


		#For an event to pass, 7 of the 9 inter comparisons must have the same sign for dpsi
		#There are 9 inter comparisons and 6 intra comparisons
		interdpsis = [KOsomaApsi - WTsomaApsi, KOsomaBpsi - WTsomaApsi, KOsomaCpsi - WTsomaApsi, 
		KOsomaApsi - WTsomaBpsi, KOsomaBpsi - WTsomaBpsi, KOsomaCpsi - WTsomaBpsi, 
		KOsomaApsi - WTsomaCpsi, KOsomaBpsi - WTsomaCpsi, KOsomaCpsi - WTsomaCpsi]
		intradpsis = [WTsomaBpsi - WTsomaApsi, WTsomaCpsi - WTsomaApsi, WTsomaCpsi - WTsomaBpsi, KOsomaBpsi - KOsomaApsi, KOsomaCpsi - KOsomaApsi, KOsomaCpsi - KOsomaBpsi]
		pos = 0
		neg = 0

		interdpsifilter = False
		for i in range(len(interdpsis)):
			if interdpsis[i] > 0:
				pos +=1
			elif interdpsis[i] < 0:
				neg +=1

		if pos >= 8:
			interdpsifilter = True


		#For an event to pass, the mean dpsi value in the inter sample comaprisons must be larger than the mean dpsi value in the intra sample comparisons
		#Doesn't add much
		intradpsifilter = False
		interdpsimean = np.mean([abs(dpsi) for dpsi in interdpsis])
		intradpsimean = np.mean([abs(dpsi) for dpsi in intradpsis])
		if interdpsimean > intradpsimean:
			intradpsifilter = True

		#To get unchanged events
		unchangedfilter = False
		unchangedcounter = 0
		for dpsi in interdpsis:
			if abs(dpsi) <= 0.01:
				unchangedcounter +=1
		if unchangedcounter >= 9 and intersig == 0:
			unchangedfilter = True


		if eventtype == 'SE':
			counter +=1
			if dpsisigfilter and dpsifilter and interdpsifilter and intradpsifilter:
				passingcounter +=1
				outfh.write(('\t').join(line) + '\n')
			if unchangedfilter:
				unchangedeventcounter +=1
				#outfh.write(('\t').join(line) + '\n')

	print counter, passingcounter, round(passingcounter / float(counter), 4) * 100, unchangedeventcounter, round(unchangedeventcounter / float(counter), 4) * 100
	infh.close()
	outfh.close()



filterdpsitable(sys.argv[1], sys.argv[2])


