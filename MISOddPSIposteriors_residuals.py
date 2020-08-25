#Created by MT on 1/21/15.

#Given a Delta psi table produced by MisoSigEventsTimecoursev2.0.py and the directory of miso outputs used to create that table
#(the directory from the --MISOdirectory flag in that script), calculate the most likely delta delta psi values for all events
#in the table as well as the pvalue for the delta delta psi having that sign.  For example, if the most likely ddPSI was -0.15,
#then the reported pvalue is the chance that the ddPSI is actually greater than or equal to 0.

#This is done using the sampled psi values in miso output files.  For each event type present in the PSI table, the chromosome
#that the event is on is first figured out by looking through all of the "chr" directories in the event type MISOdirectory subdirectories.
#The sampled psi values from two different samples (for example WTSoma_v_WTAxon) are combined to calculate a distribution of delta psis (theoretically ranging from
#-1 to 1).  The delta psi distributions from two different comparisons (e.g. WTSoma_v_WTAxon_vv_Mbnl1Soma_v_Mbnl1Axon) are then
#compared to produced ddPSI distributions (theoretically ranging from -2 to 2).


import os
import operator
import argparse
import numpy
import scipy.odr
def fitline(psitable, comparison):
    #Given a deltadeltapsicomparison (see getdeltadeltapsis, e.g. Mbnl1Soma.psi_v_Mbnl1Axon.psi_vv_WTSoma.psi_v_WTAxon.psi), 
    #go through the PSItable and get any events where BOTH samples in EITHER delta psi comparison are in the SigComparisons
    #column.  This is because the SigComparisons column could have a comparison as Mbnl1Soma.psivMbnl1Axon.psi or 
    #Mbnl1Axon.psivMbnl1Soma.psi.  Then, given those delta psi values (samplepair1 on x axis and samplepair2 on y axis, calculate
    #the slope of the line that fits those for each eventtype.

    #Delta values are always sample2 - sample 1
    #Take off '.psi'
    deltapsicomparisons = comparison.split('_vv_')
    comp1sample1 = deltapsicomparisons[0].split('_v_')[0][:-4]
    comp1sample2 = deltapsicomparisons[0].split('_v_')[1][:-4]
    comp2sample1 = deltapsicomparisons[1].split('_v_')[0][:-4]
    comp2sample2 = deltapsicomparisons[1].split('_v_')[1][:-4]
    sigdpsis = {} # {eventtype : [[comp1 delta psis from samples that pass sigcomp filter above], [comp2 delta psis from samples taht pass sigcomp filter]]}
    fits = {} # {eventtype : [slope, intercept]}
    
    #Figure out which column is the SigComparisons column and other relevant columns
    fh = open(psitable, 'r')
    header = fh.readline().strip().split('\t')
    fh.close()
    for idx, colname in enumerate(header):
        if colname == 'SigComparisons':
            sigcompidx = idx
        elif colname == comp1sample1:
            comp1sample1idx = idx
        elif colname == comp1sample2:
            comp1sample2idx = idx
        elif colname == comp2sample1:
            comp2sample1idx = idx
        elif colname == comp2sample2:
            comp2sample2idx = idx

    #Go through and get events that pass Sigcomparisons filter as described above

    fh = open(psitable, 'r')
    for line in fh:
        line = line.strip().split('\t')
        if line[0] == 'Event':
            continue
        comp1sig = False
        comp2sig = False
        eventname = line[0]
        eventtype = line[1]
        comp1sample1psi = float(line[comp1sample1idx])
        comp1sample2psi = float(line[comp1sample2idx])
        comp2sample1psi = float(line[comp2sample1idx])
        comp2sample2psi = float(line[comp2sample2idx])
        #Delta psi is always sample 2 - sample 1
        comp1deltapsi = comp1sample2psi - comp1sample1psi
        comp2deltapsi = comp2sample2psi - comp2sample1psi
        sigcomps = line[sigcompidx].split(';')
        for sigcomp in sigcomps:
            #This allows for matching to Sample1vSample2 and also Sample2vSample1
            #Can change these filters (especially deltapsi filter)
            if (comp1sample1 + 'v' + comp1sample2 in sigcomp or comp1sample2 + 'v' + comp1sample1 in sigcomp):
                comp1sig = True
            if (comp2sample1 + 'v' + comp2sample2 in sigcomp or comp2sample2 + 'v' + comp2sample1 in sigcomp):
                comp2sig = True

        #For line fitting, only consider those events that were significant in at least one of the two delta psi comparisons
        if comp1sig == True or comp2sig == True:
            if sigdpsis.has_key(eventtype) == False:
                sigdpsis[eventtype] = [[],[]]
            sigdpsis[eventtype][0].append(comp1deltapsi)
            sigdpsis[eventtype][1].append(comp2deltapsi)

    fh.close()

    #Calculate fits
    #for eventtype in sigdpsis:
        #coefficients = numpy.polyfit(sigdpsis[eventtype][0], sigdpsis[eventtype][1], 1)
        #slope = coefficients[0]
        #yint = coefficients[1]
        #fits[eventtype] = [slope, yint]
        #print 'For ddPSI comparison {0}, {1} {2} events pass BF test in at least one of the two dPSI comparisons.'.format(comparison, len(sigdpsis[eventtype][0]), eventtype) 
        #print 'The fit of the line is y = {0}x + {1}.'.format(slope, yint)

    #Alternatively, can force y intercept to be 0.
    #for eventtype in sigdpsis:
    	#x = numpy.array(sigdpsis[eventtype][0])
    	#y = numpy.array(sigdpsis[eventtype][1])
    	#x = x[:, numpy.newaxis]
    	#slope, _, _, _, = numpy.linalg.lstsq(x, y)
    	#slope = float(slope)
    	#fits[eventtype] = [slope, 0]
    	#print 'For ddPSI comparison {0}, {1} {2} events pass BF test in at least one of the two dPSI comparisons.'.format(comparison, len(sigdpsis[eventtype][0]), eventtype) 
        #print 'The fit of the line is y = {0}x.'.format(slope)

    #Alternatively, can use orthogonal distance regression to minimize x/y distance to line instead of normal regression,
    #which minimizes only y distance to line.  Seems to result in a better fit.
    for eventtype in sigdpsis:
    	x = numpy.array(sigdpsis[eventtype][0])
    	y = numpy.array(sigdpsis[eventtype][1])
    	xsd = numpy.std(x)
    	ysd = numpy.std(y)
    	def fit(B, x):
    		#This forces intercept to be 0. Can add in extra term to allow nonzero intercept: B[0]*x + B[1]
    		return B[0]*x
    	linear = scipy.odr.Model(fit)
    	data = scipy.odr.RealData(x, y, sx = xsd, sy = ysd)
    	#As a first estimate, the slope is 1
    	myodr = scipy.odr.ODR(data, linear, [1])
    	myoutput = myodr.run()
    	slope = float(myoutput.beta[0])
    	fits[eventtype] = [slope, 0]
    	print 'For ddPSI comparison {0}, {1} {2} events pass BF test in at least one of the two dPSI comparisons.'.format(comparison, len(sigdpsis[eventtype][0]), eventtype) 
        print 'The fit of the line is y = {0}x.'.format(slope)



    return fits


def getposteriors(event):
    #Given a miso file, return the frequency of each posterior psi value.
    misofh = open(event, 'r')
    posteriors = []
    posteriordist = {} #{posterior : fraction of posteriors with that value}
    for line in misofh:
        line = line.strip().split('\t')
        #skip header lines
        if line[0].startswith('#') or line[0].startswith('sampled_psi'):
            continue
        posterior = round(float(line[0].split(',')[0]), 2) #The first of the 2 given psi values.
        posteriors.append(posterior)

    misofh.close()
    if len(posteriors) == 0:
        print 'WARNING: no posterior values found for event {0}.'.format(event)
        return None

    else:
        for posterior in set(posteriors): #for each unique posterior value
            posteriordist[posterior] = posteriors.count(posterior) / float(len(posteriors))

        return posteriordist

def geteventlist(psitable):
    #Retrieve all of the eventnames and eventtypes from a PSI table produced by MISOSigEventtimecoursev2.0.py
    events = {} # {eventtype : [eventnames]}
    fh = open(psitable, 'r')
    for line in fh:
        line = line.strip().split('\t')
        if line[0] == 'Event':
            continue
        eventname = line[0]
        eventtype = line[1]
        if eventtype not in events:
            events[eventtype] = [eventname]
        elif eventtype in events:
            events[eventtype].append(eventname)

    fh.close()
    return events

def getposteriordists(events, MISOdirectory, sampledirs):
    #sampledirs are now the basenames of the sampledirs
    event_to_chrm = {} #{eventtype : {eventname : chrm}}
    sampledirs = sampledirs.split(',')

    posteriordists = {} # {eventtype : {event : {sample : {posteriorvalue : frequency}}}}
    samplenames = []
    for sampledir in sampledirs:
        if sampledir.endswith('/'):
            samplenames.append(os.path.basename(sampledir[:-1]))
        else:
            samplenames.append(os.path.basename(sampledir))

    #Figure out which chromosome directory a given event is in
    for eventtype in events:
        eventcounter = 0
        event_to_chrm[eventtype] = {}
        eventtypedir = os.path.join(os.path.abspath(MISOdirectory), eventtype)
        sample1dir = os.path.join(eventtypedir, samplenames[0])
        for eventname in events[eventtype]:
            eventcounter +=1
            if eventcounter % 100 == 0:
                print 'Finding chromosome for {0} event {1} of {2}...'.format(eventtype, eventcounter, len(events[eventtype]))
            for chrdir in os.listdir(sample1dir):
                chrdir = os.path.join(sample1dir, chrdir)
                if os.path.basename(chrdir).startswith('chr'):
                    for misofile in os.listdir(chrdir):
                        if os.path.basename(misofile) == eventname + '.miso':
                            event_to_chrm[eventtype][eventname] = os.path.basename(chrdir)

    #Get the PSI value distributions for every event in every sample
    print 'Calculating PSI value distributions...'
    for eventtype in event_to_chrm:
        eventcounter = 0
        posteriordists[eventtype] = {}
        for eventname in event_to_chrm[eventtype]:
            eventcounter +=1
            if eventcounter % 100 == 0:
                print 'Calculating PSI value for {0} event {1}.'.format(eventtype, eventcounter)
            #First, populate posteriordists dictionary
            posteriordists[eventtype][eventname] = {}
            for samplename in samplenames:
                posteriordists[eventtype][eventname][samplename] = {}
            chrm = event_to_chrm[eventtype][eventname]
            misofiles = [os.path.join(os.path.abspath(MISOdirectory), eventtype, sampledir, chrm, eventname + '.miso') for sampledir in sampledirs]
            for misofile in misofiles:
                misofile = os.path.abspath(misofile)
                samplename = misofile.split('/')[-3]
                posteriordist = getposteriors(misofile)
                posteriordists[eventtype][eventname][samplename] = posteriordist

   
    return posteriordists


def getdeltapsidists(posteriordists, deltapsicomparisons):
    #Given a dictionary of posterior distributions, calculate the deltapsi distributions between sample1 and sample2
    #This works best if sample1 is a soma sample and sample2 is an axon sample as deltapsi will be calculated as sample2 - sample1
    #Sample names should be the basenames of the MISO output directory, e.g. 'Mbnl1Soma.psi', so that they will match the keys
    #in the posteriordist dictionary returned by getposteriordists
    #Comparisons are a comma separated list of pairwise comparisons to be made with '_v_' between individual samples in a comparison, e.g.:
    #Mbnl1Soma.psi_v_Mbnl1Axon.psi,Mbnl2Soma.psi_v_Mbnl2Axon.psi

    #posteriordists = {eventtype : {event : {sample : {posteriorvalue : frequency}}}}

    comparisons = deltapsicomparisons.split(',') #Mbnl1Soma.psivMbnl1Axon.psi,Mbnl2Soma.psivMbnl2Axon.psi becomes ['Mbnl1Soma.psi_v_Mbnl1Axon.psi','Mbnl2Soma.psi_v_Mbnl2Axon.psi']
    deltapsidists = {} # {eventtype: {comparison : {event : {deltapsi : P}}}}

    print 'Calculating delta PSI distributions...'

    for eventtype in posteriordists:
        print 'Calculating delta PSI distributions for {0} events.'.format(eventtype)
        deltapsidists[eventtype] = {}
        for comparison in comparisons:
            deltapsidists[eventtype][comparison] = {}
            sample1 = comparison.split('_v_')[0]
            sample2 = comparison.split('_v_')[1]
            for event in posteriordists[eventtype]:
                deltapsidists[eventtype][comparison][event] = {}
                for posteriorvalue_s2 in posteriordists[eventtype][event][sample2]:
                    for posteriorvalue_s1 in posteriordists[eventtype][event][sample1]:
                        #For every possible sample1 posteriorvalue/ sample2 posteriorvalue combination, calculate deltapsi
                        #The P of that deltapsi is the product of the frequencies of the psi values in their respective samples
                        #If that deltapsi hasn't been seen before for this event, that's now the P for it
                        #If it has been seen before for this event, add it to the preexisting P
                        deltapsi = round(posteriorvalue_s2 - posteriorvalue_s1, 2)
                        s2freq = posteriordists[eventtype][event][sample2][posteriorvalue_s2]
                        s1freq = posteriordists[eventtype][event][sample1][posteriorvalue_s1]
                        P = s2freq * s1freq
                        if deltapsi not in deltapsidists[eventtype][comparison][event]:
                            deltapsidists[eventtype][comparison][event][deltapsi] = P
                        elif deltapsi in deltapsidists[eventtype][comparison][event]:
                            newP = deltapsidists[eventtype][comparison][event][deltapsi] + P
                            deltapsidists[eventtype][comparison][event][deltapsi] = newP

    return deltapsidists


def getdeltadeltapsidists(deltapsidists, deltadeltapsicomparisons, psitable):
    #Given a dictionary of deltapsi distributions, calculate the deltadeltapsi distributions between two delta psi distributions (deltapsidist1 and deltapsidist2).
    #This works best if deltapsidist2 is "treatment" and deltapsidist1 is "control" as deltadeltapsi will be calculated as deltapsi2 - deltapsi1
    #Deltadeltapsi comparison names are deltapsi comparison names separated by "_vv_"
    #E.g. Mbnl1Soma.psi_v_Mbnl1Axon.psi_vv_WTSoma.psi_v_WTAxon.psi

    #Later (05/22/15) added the ability to consider "residual" delta delta psis.  This is because often there is a uniform increase/decrease in delta psis between
    #two samples.  To try to correct for that, fit a line to a deltapsi/deltapsi scatter (between two conditions).  Then, given a delta psi in one condition and the
    #equation of that line, you can calculate the expected delta psi in the other condition.  Any deviation from that is the "residual", which will be treated as 
    #a delta delta psi.  Since every residual value is directly related to one and only one delta psi value pair, the probability of getting that residual
    #is the same probability as getting a "naive" delta delta psi, i.e. just subtracting the two delta psis.  This means that the probabilities for residuals can be
    #calculated just as the probabilities for delta psis were calculated.

    
    comparisons = deltadeltapsicomparisons.split(',') 
    #Mbnl1Soma.psi_v_Mbnl1Axon.psi_vv_WTSoma.psi_v_WTAxon.psi,Mbnl2Soma.psi_v_Mbnl2Axon.psi_vv_WTSoma.psi_v_WTAxon.psi becomes 
    #['Mbnl1Soma.psi_v_Mbnl1Axon.psi_vv_WTSoma.psi_v_WTAxon.psi' , 'Mbnl2Soma.psi_v_Mbnl2Axon.psi_vv_WTSoma.psi_v_WTAxon.psi']

    #deltapsidists = {eventtype: {comparison : {event : {deltapsi : P}}}}
    #deltadeltapsidists = {} # {eventtype: {deltadeltapsicomparison : {event : {deltadeltapsi : P}}}}
    residualdists = {} # {eventtype : {deltadeltapsicomparison : {event : {residual : P}}}}

    for eventtype in deltapsidists:
        #deltadeltapsidists[eventtype] = {}
        residualdists[eventtype] = {}
        for comparison in comparisons:
            fits = fitline(psitable, comparison)
            slope = fits[eventtype][0]
            yint = fits[eventtype][1]
            #deltadeltapsidists[eventtype][comparison] = {}
            residualdists[eventtype][comparison] = {}
            samplepair1 = comparison.split('_vv_')[0]
            samplepair2 = comparison.split('_vv_')[1]
            print 'Performing {0} vs {1} delta delta psi calculations.'.format(samplepair1, samplepair2)

            samplepair1events = deltapsidists[eventtype][samplepair1]
            samplepair2events = deltapsidists[eventtype][samplepair2]

            for event in samplepair1events:
                #deltadeltapsidists[eventtype][comparison][event] = {}
                residualdists[eventtype][comparison][event] = {}
                samplepair1deltapsis = samplepair1events[event]
                samplepair2deltapsis = samplepair2events[event]
                for samplepair1deltapsi in samplepair1deltapsis: #compare all possible pairwise combinations between deltapsis of the two sample pairs
                    for samplepair2deltapsi in samplepair2deltapsis:
                        expecteddeltapsi = (samplepair1deltapsi * slope) + yint
                        residual = samplepair2deltapsi - expecteddeltapsi
                        deltadeltapsi = round(samplepair2deltapsi - samplepair1deltapsi, 2)
                        samplepair1freq = samplepair1deltapsis[samplepair1deltapsi]
                        samplepair2freq = samplepair2deltapsis[samplepair2deltapsi]
                        P = samplepair1freq * samplepair2freq #prob of getting that ddPSI, therefore also prob of having that "residual"
                        #print event, residual, P
                        #if deltadeltapsi not in deltadeltapsidists[eventtype][comparison][event]: #if we haven't seen this deltadeltapsi for this event in this comparison before
                            #deltadeltapsidists[eventtype][comparison][event][deltadeltapsi] = P
                        #elif deltadeltapsi in deltadeltapsidists[eventtype][comparison][event]: #if we have seen this deltadeltapsi for this event in this comparison before
                            #newP = deltadeltapsidists[eventtype][comparison][event][deltadeltapsi] + P
                            #deltadeltapsidists[eventtype][comparison][event][deltadeltapsi] = newP

                        if residual not in residualdists[eventtype][comparison][event]:
                            residualdists[eventtype][comparison][event][residual] = P
                        elif residual in residualdists[eventtype][comparison][event]:
                            newP = residualdists[eventtype][comparisons][event][residual] + P
                            residualdists[eventtype][comparison][event][residual] = newP

    return residualdists
    #return deltadeltapsidists

def summarizedeltadeltapsidists(deltadeltapsidists):
    #deltadeltapsidists = {eventtype: {deltadeltapsicomparison : {event : {deltadeltapsi : P}}}}

    summarizeddists = {} #{eventtype : {event : {comparison : {mostlikelydeltadeltapsi : p value on the ddPSI having sign of mostlikely ddPSI}}}}
    print 'Summarizing delta delta PSI distributions...'

    #initialize all dictionaries in summarizeddists
    for eventtype in deltadeltapsidists:
        events = []
        for deltapsicomparison in deltadeltapsidists[eventtype]:
            for event in deltadeltapsidists[eventtype][deltapsicomparison]:
                if event not in events:
                    events.append(event)
        print len(events)
        comparisons = []
        for comparison in deltadeltapsidists[eventtype]:
            if comparison not in comparisons:
                comparisons.append(comparison)
        print len(comparisons)

        summarizeddists[eventtype] = {}


    
        for event in events:
            summarizeddists[eventtype][event] = {}
            for comparison in comparisons:
                summarizeddists[eventtype][event][comparison] = {}
                ddPSIfreqs = deltadeltapsidists[eventtype][comparison][event]
                mostlikelyddpsi = max(ddPSIfreqs.iteritems(), key = operator.itemgetter(1))[0] #Get the delta delta psi with the largest frequency.
                if mostlikelyddpsi == 0:  #If ddPSI = 0, then the pvalue is 1
                    pvalue = 1.0
                elif mostlikelyddpsi > 0:  #If ddPSI is positive, the pvalue is the summed frequency of all ddPSIs less than or equal to 0
                    pvaluefreqs = []
                    for ddPSI in ddPSIfreqs:
                        if ddPSI <= 0:
                            pvaluefreqs.append(ddPSIfreqs[ddPSI])
                        
                elif mostlikelyddpsi < 0: #If ddPSI is negative, the pvalue is the summed frequence of all ddPSIs greater than or equal to 0
                    pvaluefreqs = []
                    for ddPSI in ddPSIfreqs:
                        if ddPSI >= 0:
                            pvaluefreqs.append(ddPSIfreqs[ddPSI])
                
                if pvaluefreqs:
                    pvalue = sum(pvaluefreqs)
                elif len(pvaluefreqs) == 0:
                    pvalue = 0.0000000001
                mostlikelyddpsi = round(mostlikelyddpsi, 4)

                summarizeddists[eventtype][event][comparison][mostlikelyddpsi] = pvalue

    return summarizeddists


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PSItable', required = True, type = str, help = 'PSI table produced by MISOSigEventsTimecourse.py')
    parser.add_argument('--MISOdirectory', required = True, type = str, help = 'The top-level directory of MISO results. Contains subdirectories for each eventtype, which themselves contain directories for each sample. Identical to MISOdirectory flag in MISOSigEventsTimecourse.py')
    parser.add_argument('--sampledirs', required = True, type = str, help = 'Comma separated list of sample directory basenames. For example, WTSoma.psi,WTAxon.psi,Mbnl1Soma.psi,Mbnl1Axon.psi. Each of these directories contains subdirectories for chromosomes.')
    parser.add_argument('--deltapsicomparisons', required = True, type = str, help = 'Comma separated list of sample comparisons for delta psi calculations. In each comparison, each sample is separated by \'_v_\', e.g. Mbnl1Soma.psi_v_Mbnl1Axon.psi,Mbnl2Soma.psi_v_Mbnl2Axon.psi. Delta psis are calculated as sample2 - sample1.')
    parser.add_argument('--deltadeltapsicomparisons', required = True, type = str,
                        help = 'Comma separated list of sample pairs for delta delta psi calculations.  In each comparison, each sample pair is separated by \'_vv_\', e.g. WTSoma.psi_v_WTAxon.psi_vv_Mbnl1Soma.psi_v_Mbnl1Axon.psi,WTSoma.psi_v_WTAxon.psi_vv_Mbnl2Soma.psi_v_Mbnl2Axon.psi. Sample pairs must be in same order as they are presented in deltapsicomparisons.  Delta delta psis are calculated as samplepair2 - samplepair1.')
    parser.add_argument('--outfile', type = str, required = True, help = 'Output file. Similar to starting PSI table but with extra fields.')
    args = parser.parse_args()

events = geteventlist(args.PSItable)
posteriordists = getposteriordists(events, args.MISOdirectory, args.sampledirs)
deltapsidists = getdeltapsidists(posteriordists, args.deltapsicomparisons)
deltadeltapsidists = getdeltadeltapsidists(deltapsidists, args.deltadeltapsicomparisons, args.PSItable)             
summarizeddists = summarizedeltadeltapsidists(deltadeltapsidists)
dPSIcomps = args.deltapsicomparisons.split(',')
ddPSIcomps = args.deltadeltapsicomparisons.split(',')
numddPSIcomps = len(ddPSIcomps)


#summarizeddists = {eventtype : {event : {comparison : {mostlikelydeltadeltapsi : p value on the ddPSI having sign of mostlikely ddPSI}}}}

infh = open(args.PSItable, 'r')
outfh = open(args.outfile, 'w')


for line in infh:
    line = line.strip().split('\t')
    event = line[0]
    eventtype = line[1]
    if line[0] == 'Event': #if this is the header
        outfh.write(('\t').join(line) + '\t' + ('\t').join(ddPSIcomps))
        for i in range(numddPSIcomps):
            outfh.write('\t' + str(ddPSIcomps[i]) + '_pvalue')
        outfh.write('\n')
    else:
        outfh.write(('\t').join(line))
        for i in range(numddPSIcomps):
            outfh.write('\t' + str(summarizeddists[eventtype][event][ddPSIcomps[i]].keys()[0]))
        for i in range(numddPSIcomps):
            outfh.write('\t' + str(summarizeddists[eventtype][event][ddPSIcomps[i]].values()[0]))
        outfh.write('\n')

infh.close()
outfh.close()
