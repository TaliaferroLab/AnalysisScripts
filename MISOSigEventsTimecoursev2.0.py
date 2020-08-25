#This script takes a MISO output directory of a particular organization (described below), and looks
#through comparisons between samples in that directory for events that are significant in at least one
#comparison.  It then takes those events and gets PSI values for them in each sample.  Psi values
#are only reported if they were calculated for every sample.

#Version 2.0 also outputs the comparisons in which each event was significant.  Version 1.0 does not do this.

#Usage: python MISOSigEventsTimecourse.py --help

import sys
import os
import argparse
import json
from collections import defaultdict

#Read count filters.  Change if necessary.
readcountfilters = {} #{eventtype : [inclusion_counts (1,0), exclusion_counts (0,1), total_counts ((1,0) + (0,1))]}
readcountfilters['AFE'] = [10, 10, 25]
readcountfilters['ALE'] = [20, 20, 50]
readcountfilters['MXE'] = [1, 1, 10]
readcountfilters['RI'] = [1, 1, 10]
readcountfilters['SE'] = [1, 1, 10]
readcountfilters['TandemUTR'] = [10, 0, 50]
#readcountfilters['SE_NMDINC'] = [1, 1, 10]
#readcountfilters['SE_NMDEXC'] = [1, 1, 10]    


#This function makes a nested dictionary
def makehash():
    return defaultdict(makehash)

#This function gets all events that are significant in at least one comparison
def MISOcompare(MISOdirectory): 
    #Main MISO Output directory. Contains subdirectories for all event types.  These
    #subdirectories contains directories called "comparisons" which themselves contain
    #directories for each pairwise comparison (e.g. CADAxon1.psi_vs_CADAxon2.psi)
    #which themselves contain a directory called bayes-factors which has in it a file 
    #called <comparison>.miso_bf.
    eventtypes = ['AFE','ALE','MXE','RI','SE','TandemUTR'] #change this according to event type you are interested in
    eventtypedirs = [] #list of all <eventtype>/comparisons directories in MISOdirectory
    outlists = []
    sigevents = {} # {eventtype : [list of sig events in at least one comparison]}

    #Populate sigevents
    for eventtype in eventtypes:
        sigevents[eventtype] = {}

    #Get eventtype dirs
    for eventtype in eventtypes:
        eventtypedir = os.path.join(os.path.abspath(MISOdirectory), eventtype, 'comparisons')
        if os.path.exists(eventtypedir):
            eventtypedirs.append(eventtypedir)

    for eventtypedir in eventtypedirs:
        #eventtypedir = <MISOdirectory>/<eventtype>/comparisons
        #One directory up would be the <eventtype> directory
        eventtype = os.path.basename(os.path.abspath(os.path.join(eventtypedir, '..')))
        #eventtype = 'AFE', 'ALE', etc.
        comparisons = []
        #Get list of comparisons to consider
        for directory, dirnames, filenames in os.walk(eventtypedir):
            for dirname in dirnames:
                if '_vs_' in dirname:
                    comparisons.append(dirname)

        for comparison in comparisons:
            comparison_results = [] #[sample1, sample2, eventtype, #events, percent events that are sig.]
            eventcounter = 0
            sigeventcounter = 0
            sample1 = comparison.split('.')[0]
            #sample2 = comparison.split('.')[1].split('_')[2]
            sample2 = comparison.split('_vs_')[1].split('.psi')[0]
            
            #Get bayes factor file
            bffile = os.path.join(os.path.abspath(eventtypedir), comparison, 'bayes-factors', comparison + '.miso_bf')
            bffh = open(bffile, 'r')
            for line in bffh:
                line = line.strip().split('\t')
                if line[0] == 'event_name': #skip header
                    continue
                eventname = line[0]
                numberofisoforms = len(line[1].split(',')) + 1
                sample1countspass = False
                sample2countspass = False #reset these values
                #Only consider events with 2 isoforms
                if numberofisoforms == 2:
                    eventcounter +=1
                    inclusionindex = None
                    exclusionindex = None #reset these values
                    deltapsi = float(line[7])
                    bayesfactor = float(line[8])
                    sample1counts = line[10]
                    sample1counts = sample1counts.replace(',(', ';(') #replace comma between classes with semicolon
                    sample1countsd = dict(item.split(':') for item in sample1counts.split(';'))
                    sample2counts = line[12]
                    sample2counts = sample2counts.replace(',(', ';(')
                    sample2countsd = dict(item.split(':') for item in sample2counts.split(';'))
                    for countclass in sample1countsd:
                        sample1countsd[countclass] = int(sample1countsd[countclass]) #change read counts to integers
                    for countclass in sample2countsd:
                        sample2countsd[countclass] = int(sample2countsd[countclass])
                    if '(1,0)' not in sample1countsd:
                        sample1countsd['(1,0)'] = 0
                    if '(0,1)' not in sample1countsd:
                        sample1countsd['(0,1)'] = 0
                    if '(1,0)' not in sample2countsd:
                        sample2countsd['(1,0)'] = 0
                    if '(0,1)' not in sample2countsd:
                        sample2countsd['(0,1)'] = 0
                    combinedsample1counts = sample1countsd['(1,0)'] + sample1countsd['(0,1)']
                    combinedsample2counts = sample2countsd['(1,0)'] + sample2countsd['(0,1)']
                    readfilters = readcountfilters[eventtype] #get read coverage filters
                    
                    if (sample1countsd['(1,0)'] >= readfilters[0] and sample1countsd['(0,1)'] >= readfilters[1] 
                        and combinedsample1counts >= readfilters[2]):
                        sample1countspass = True

                    if (sample2countsd['(1,0)'] >= readfilters[0] and sample2countsd['(0,1)'] >= readfilters[1] 
                        and combinedsample2counts >= readfilters[2]):
                        sample2countspass = True

                    #Does it pass all filters?
                    if abs(deltapsi) >= 0 and bayesfactor >= 10 and sample1countspass and sample2countspass:
                        sigeventcounter +=1
                        if eventname not in sigevents[eventtype]:
                            sigevents[eventtype][eventname] = ['{0}v{1}'.format(sample1, sample2)]
                        elif eventname in sigevents[eventtype]:
                            sigevents[eventtype][eventname].append('{0}v{1}'.format(sample1, sample2))

            print '{0} significant events in {1} {2}.'.format(sigeventcounter, eventtype, comparison)

    #Remove duplicates in each list of significant events
    for eventtype in sigevents:
        print 'In total, {0} {1} events are significant it at least one comparison.'.format(len(sigevents[eventtype]), eventtype)

    return sigevents

def getPSIs(MISOdirectory, sigevents):
    #MISOdirectory again is the top level directory. Contains subdirectories of each eventtype containing MISO output.
    #Get PSI values in each sample for those events that you retrieved using MISOcompare
    #These should all be 2 isoform events, so getting PSI value should be straightforward
    eventtypes = ['AFE','ALE','MXE','RI','SE','TandemUTR'] #change this according to event type you are interested in
    samples = []
    sigcomparisondict = {} #{eventtype : {event : [list of sig comparisons]}}
    PSIs = makehash() # {eventtype : {event : {sample : psivalue}}}

    for eventtype in eventtypes:
        sigcomparisondict[eventtype] = {}
        sampledirs = [] #e.g. CADAxon1.psi
        eventtypedir = os.path.join(os.path.abspath(MISOdirectory), eventtype)
        for sampledir in os.listdir(eventtypedir):
            if sampledir.endswith('.psi'):
                sampledirs.append(os.path.join(os.path.abspath(MISOdirectory), eventtype, sampledir))
                samples.append(os.path.basename(sampledir).split('.')[0])

        for sampledir in sampledirs:
            sample = os.path.basename(sampledir).split('.')[0]
            print 'Retrieving PSI values for {0} {1}...'.format(sample, eventtype)
            misosummaryfile = os.path.join(os.path.abspath(MISOdirectory), eventtype, sampledir, 'summary_output', 'summary', '{0}.miso_summary'.format(os.path.basename(sampledir)))
            misosummaryfh = open(misosummaryfile, 'r')
            for line in misosummaryfh:
                line = line.strip().split('\t')
                eventname = line[0]
                if eventname == 'event_name': #skip header
                    continue
                if eventname in sigevents[eventtype]:
                    PSI = float(line[1])
                    sigcomparisons = sigevents[eventtype][eventname]
                    PSIs[eventtype][eventname][sample] = PSI
                    sigcomparisondict[eventtype][eventname] = sigcomparisons

    samples = list(set(samples))

    return PSIs, samples, sigcomparisondict
                
                
            
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--MISOdirectory', type = str, help = 'Main MISO Output directory. Contains subdirectories for all event types.  These subdirectories contains directories called "comparisons" which themselves contain directories for each pairwise comparison (e.g. CADAxon1.psi_vs_CADAxon2.psi) which themselves contain a directory called bayes-factors. Required.')
    parser.add_argument('--outfile', type = str, help = 'Output file. Required.')
    args = parser.parse_args()
    sigevents = MISOcompare(args.MISOdirectory)
    
    PSIs, samples, sigcomparisondict = getPSIs(args.MISOdirectory, sigevents)
    outfh = open(args.outfile , 'w')
    outfh.write('Event' + '\t' + 'EventType' + '\t' + ('\t').join(sorted(samples)) + '\t' 'SigComparisons' + '\n')
    for eventtype in PSIs:
        for event in PSIs[eventtype]:
            if len(PSIs[eventtype][event]) == len(samples): #only consider events with a PSI value for every sample
                outfh.write(event + '\t' + eventtype)
                for sample in sorted(PSIs[eventtype][event]): #this is important so that PSI values are written in correct order
                    #The length of sigcomparisondict[eventtype][event] is how many comparisons this event was significant in
                    outfh.write('\t' + str(PSIs[eventtype][event][sample]))
                outfh.write('\t' + (';').join(sigcomparisondict[eventtype][event]))
                outfh.write('\n')
            
            

    outfh.close()
