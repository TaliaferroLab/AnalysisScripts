#Takes a harsh view of whether or not an event "passes" and is significant.  Assumes there are three replicates
#of two cell lines(N2A and CAD).  There are the three "test" comparisons (soma1 vs axon1, soma2 vs axon2, soma3 vs axon3)
#and six "ctrl" comparisons (soma 1 vs soma2, soma1 vs soma3, soma2 vs soma3, axon1 vs axon2, axon2 vs axon3, axon1 vs axon3).
#If an event passes filters in a comparison with a positive delta psi, it is designated "positive."
#If an event passes filters in a comparison with a negative delta psi, it is designated "negative."
#If an event doesn't pass filters, it is designated "fail."

#An event passes if it is "positive" in at least 2 out of 3 or "negative" in 2 out of 3 "test" comparisons
#or 2 out of 3 randomly chosen (from the 6) control comparisons.

#If --celllines == together, then an event must pass in both cell lines (with the same sign) to pass.

import sys
import os
import argparse
import random


#Read count filters.  Change if necessary.
readcountfilters = {} #{eventtype : [inclusion_counts (1,0), exclusion_counts (0,1), total_counts ((1,0) + (0,1))]}
readcountfilters['AFE'] = [50, 50, 100]
readcountfilters['ALE'] = [50, 50, 100]
readcountfilters['MXE'] = [1, 1, 10]
readcountfilters['RI'] = [1, 1, 10]
readcountfilters['SE'] = [10, 10, 25]
readcountfilters['TandemUTR'] = [50, 0, 100]    

def MISOcompare(MISOdirectory): 
    #Main MISO Output directory. Contains subdirectories for all event types.  These
    #subdirectories contains directories called "comparisons" which themselves contain
    #directories for each pairwise comparison (e.g. CADAxon1.psi_vs_CADAxon2.psi)
    #which themselves contain a directory called bayes-factors.
    eventtypes = ['AFE', 'ALE', 'MXE', 'RI', 'SE', 'TandemUTR']
    eventtypedirs = [] #list of all <eventtype>/comparisons directories in MISOdirectory
    outlists = []
    eventdict = {} # {eventtype : {event : {'N2Atest' : [N2Atestresults (positive/negative/fail)], 'N2Actrl' : [N2Actrlresults], 'CADtest' : [CADtestresults], 'CADctrl' : [CADctrlresults]}}}
    #Get eventtype dirs
    for eventtype in eventtypes:
        eventtypedir = os.path.join(os.path.abspath(MISOdirectory), eventtype, 'comparisons')
        if os.path.exists(eventtypedir):
            eventtypedirs.append(eventtypedir)

    for eventtypedir in eventtypedirs:
        #eventtypedir = <MISOdirectory>/<eventtype>/comparisons
        #One directory up would be the <eventtype> directory
        eventtype = os.path.basename(os.path.abspath(os.path.join(eventtypedir, '..')))
        eventdict[eventtype] = {}
        #eventtype = 'AFE', 'ALE', etc.
        print 'Going through {0} events...'.format(eventtype)
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
            sample2 = comparison.split('.')[1].split('_')[2]
            comparisontype = ''
            #Explicitly say which comparisons you want
            if 'N2ASoma1' in sample1 and 'N2AAxon1' in sample2 and sample1 != sample2:
                comparisontype = 'N2Atest'
            elif 'N2ASoma2' in sample1 and 'N2AAxon2' in sample2 and sample1 != sample2:
                comparisontype = 'N2Atest'
            elif 'N2ASoma3' in sample1 and 'N2AAxon3' in sample2 and sample1 != sample2:
                comparisontype = 'N2Atest'
            elif 'CADSoma1' in sample1 and 'CADAxon1' in sample2 and sample1 != sample2:
                comparisontype = 'CADtest'
            elif 'CADSoma2' in sample1 and 'CADAxon2' in sample2 and sample1 != sample2:
                comparisontype = 'CADtest'
            elif 'CADSoma3' in sample1 and 'CADAxon3' in sample2 and sample1 != sample2:
                comparisontype = 'CADtest'
            elif 'N2ASoma1' in sample1 and 'N2ASoma2' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'N2ASoma1' in sample1 and 'N2ASoma3' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'N2ASoma2' in sample1 and 'N2ASoma3' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'N2AAxon1' in sample1 and 'N2AAxon2' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'N2AAxon1' in sample1 and 'N2AAxon3' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'N2AAxon2' in sample1 and 'N2AAxon3' in sample2 and sample1 != sample2: 
                comparisontype = 'N2Actrl'
            elif 'CADSoma1' in sample1 and 'CADSoma2' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            elif 'CADSoma1' in sample1 and 'CADSoma3' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            elif 'CADSoma2' in sample1 and 'CADSoma3' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            elif 'CADAxon1' in sample1 and 'CADAxon2' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            elif 'CADAxon1' in sample1 and 'CADAxon3' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            elif 'CADAxon2' in sample1 and 'CADAxon3' in sample2 and sample1 != sample2:
                comparisontype = 'CADctrl'
            else: #if its not one of these specific comparisons, skip it
                continue

    #Get bayes factor file
            bffile = os.path.join(os.path.abspath(eventtypedir), comparison, 'bayes-factors', comparison + '.miso_bf')
            bffh = open(bffile, 'r')
            for line in bffh:
                line = line.strip().split('\t')
                if line[0] == 'event_name': #skip header
                    continue
                eventname = line[0]
                if eventdict[eventtype].has_key(eventname) == False:
                    eventdict[eventtype][eventname] = {}
                if eventdict[eventtype][eventname].has_key(comparisontype) == False:
                    eventdict[eventtype][eventname][comparisontype] = []
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
                    if deltapsi >= 0.1 and bayesfactor >= 10 and sample1countspass and sample2countspass:
                        sigeventcounter +=1
                        eventdict[eventtype][eventname][comparisontype].append('positive')
                    elif deltapsi <= -0.1 and bayesfactor >= 10 and sample1countspass and sample2countspass:
                        sigeventcounter +=1
                        eventdict[eventtype][eventname][comparisontype].append('negative')
                    else:
                        eventdict[eventtype][eventname][comparisontype].append('fail')

            bffh.close()

    return eventdict

def summarizeeventdict(eventdict):
    #eventdict = {} # {eventtype : {event : {'N2Atest' : [N2Atestresults (positive/negative/fail)], 'N2Actrl' : [N2Actrlresults], 'CADtest' : [CADtestresults], 'CADctrl' : [CADctrlresults]}}}
    summarized = {} # {eventtype : {'N2Atest' : [passed, failed]}, {'N2Actrl' : [passed, failed]}, {'CADtest' : [passed, failed]}, {'CADctrl' : [passed, failed]}}
    for eventtype in eventdict:
        print 'Summarizing {0} results...'.format(eventtype)
        summarized[eventtype] = {}
        summarized[eventtype]['N2Atest'] = [0, 0]
        summarized[eventtype]['N2Actrl'] = [0, 0]
        summarized[eventtype]['CADtest'] = [0, 0]
        summarized[eventtype]['CADctrl'] = [0, 0]
        for event in eventdict[eventtype]:
            #Only look at events where a comparison was made between all possible samples
            #This exculdes events that were not expressed in all samples
            if eventdict[eventtype][event].has_key('N2Atest') == False or eventdict[eventtype][event].has_key('N2Actrl') == False or eventdict[eventtype][event].has_key('CADtest') == False or eventdict[eventtype][event].has_key('CADctrl') == False:
                continue
            if len(eventdict[eventtype][event]['N2Atest']) == 3 and len(eventdict[eventtype][event]['N2Actrl']) == 6 and len(eventdict[eventtype][event]['CADtest']) == 3 and len(eventdict[eventtype][event]['CADctrl']) == 6:
                N2Atestresults = eventdict[eventtype][event]['N2Atest']
                N2Actrlresults = eventdict[eventtype][event]['N2Actrl']
                CADtestresults = eventdict[eventtype][event]['CADtest']
                CADctrlresults = eventdict[eventtype][event]['CADctrl']
                if N2Atestresults.count('positive') >=2 or N2Atestresults.count('negative') >=2:
                    summarized[eventtype]['N2Atest'][0] +=1
                else:
                    summarized[eventtype]['N2Atest'][1] +=1
                randomN2Actrls = random.sample(N2Actrlresults, 3)
                if randomN2Actrls.count('positive') >=2 or randomN2Actrls.count('negative') >=2:
                    summarized[eventtype]['N2Actrl'][0] +=1
                else:
                    summarized[eventtype]['N2Actrl'][1] +=1
                if CADtestresults.count('positive') >=2 or CADtestresults.count('negative') >=2:
                    summarized[eventtype]['CADtest'][0] +=1
                else:
                    summarized[eventtype]['CADtest'][1] +=1
                randomCADctrls = random.sample(CADctrlresults, 3)
                if randomCADctrls.count('positive') >=2 or randomCADctrls.count('negative') >=2:
                    summarized[eventtype]['CADctrl'][0] +=1
                else:
                    summarized[eventtype]['CADctrl'][1] +=1

    return summarized

def summarizeeventdict_bothcelllines(eventdict):
    #This is like the above function, but in order to pass, an event must have passed filters in BOTH cell lines
    #eventdict = {} # {eventtype : {event : {'N2Atest' : [N2Atestresults (positive/negative/fail)], 'N2Actrl' : [N2Actrlresults], 'CADtest' : [CADtestresults], 'CADctrl' : [CADctrlresults]}}}
    summarized = {} # {eventtype : {'test' : [passed, failed]}, {'ctrl' : [passed, failed]}}
    for eventtype in eventdict:
        print 'Summarizing {0} results...'.format(eventtype)
        summarized[eventtype] = {}
        summarized[eventtype]['test'] = [0, 0]
        summarized[eventtype]['ctrl'] = [0, 0]
        for event in eventdict[eventtype]:
            #Only look at events where a comparison was made between all possible samples
            #This exculdes events that were not expressed in all samples
            if eventdict[eventtype][event].has_key('N2Atest') == False or eventdict[eventtype][event].has_key('N2Actrl') == False or eventdict[eventtype][event].has_key('CADtest') == False or eventdict[eventtype][event].has_key('CADctrl') == False:
                continue
            if len(eventdict[eventtype][event]['N2Atest']) == 3 and len(eventdict[eventtype][event]['N2Actrl']) == 6 and len(eventdict[eventtype][event]['CADtest']) == 3 and len(eventdict[eventtype][event]['CADctrl']) == 6:
                N2Atestresults = eventdict[eventtype][event]['N2Atest']
                N2Actrlresults = eventdict[eventtype][event]['N2Actrl']
                CADtestresults = eventdict[eventtype][event]['CADtest']
                CADctrlresults = eventdict[eventtype][event]['CADctrl']
                if (N2Atestresults.count('positive') >=2 and CADtestresults.count('positive') >=2) or (N2Atestresults.count('negative') >=2 and CADtestresults.count('negative') >=2):
                    summarized[eventtype]['test'][0] +=1
                else:
                    summarized[eventtype]['test'][1] +=1
                randomN2Actrls = random.sample(N2Actrlresults, 3)
                randomCADctrls = random.sample(CADctrlresults, 3)
                if (randomN2Actrls.count('positive') >=2 and randomCADctrls.count('positive') >=2) or (randomN2Actrls.count('negative') >=2 and randomCADctrls.count('negative') >=2):
                    summarized[eventtype]['ctrl'][0] +=1
                else:
                    summarized[eventtype]['ctrl'][1] +=1

    return summarized


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--MISOdirectory', type = str, help = 'Main MISO Output directory. Contains subdirectories for all event types.  These subdirectories contains directories called "comparisons" which themselves contain directories for each pairwise comparison (e.g. CADAxon1.psi_vs_CADAxon2.psi) which themselves contain a directory called bayes-factors. Required.')
    parser.add_argument('--celllines', type = str, choices = ['separate', 'together'], help = 'Should cell lines be considered separately or together? If together, an event must pass filters in BOTH cell lines to pass.')
    parser.add_argument('--outfile', type = str, help = 'Output file. Required.')
    args = parser.parse_args()

    eventdict = MISOcompare(args.MISOdirectory)

    if args.celllines == 'separate':
        summarized = summarizeeventdict(eventdict)
        outfh = open(args.outfile, 'w')
        outfh.write(('\t').join(['Eventtype','Comparisontype','Number_passed','Number_failed','Fraction_passed']) + '\n')
        for eventtype in summarized:
            N2Atestpassed = summarized[eventtype]['N2Atest'][0]
            N2Atestfailed = summarized[eventtype]['N2Atest'][1]
            N2Actrlpassed = summarized[eventtype]['N2Actrl'][0]
            N2Actrlfailed = summarized[eventtype]['N2Actrl'][1]
            CADtestpassed = summarized[eventtype]['CADtest'][0]
            CADtestfailed = summarized[eventtype]['CADtest'][1]
            CADctrlpassed = summarized[eventtype]['CADctrl'][0]
            CADctrlfailed = summarized[eventtype]['CADctrl'][1]
            N2Atestfracpassed = N2Atestpassed / float(N2Atestpassed + N2Atestfailed)
            N2Actrlfracpassed = N2Actrlpassed / float(N2Actrlpassed + N2Actrlfailed)
            CADtestfracpassed = CADtestpassed / float(CADtestpassed + CADtestfailed)
            CADctrlfracpassed = CADctrlpassed / float(CADctrlpassed + CADctrlfailed)
            outfh.write(('\t').join([eventtype, 'N2Atest', str(N2Atestpassed), str(N2Atestfailed), str(N2Atestfracpassed)]) + '\n')
            outfh.write(('\t').join([eventtype, 'N2Actrl', str(N2Actrlpassed), str(N2Actrlfailed), str(N2Actrlfracpassed)]) + '\n')
            outfh.write(('\t').join([eventtype, 'CADtest', str(CADtestpassed), str(CADtestfailed), str(CADtestfracpassed)]) + '\n')
            outfh.write(('\t').join([eventtype, 'CADctrl', str(CADctrlpassed), str(CADctrlfailed), str(CADctrlfracpassed)]) + '\n')
        outfh.close()

    elif args.celllines == 'together':
        summarized = summarizeeventdict_bothcelllines(eventdict)
        outfh = open(args.outfile, 'w')
        outfh.write(('\t').join(['Eventtype','Comparisontype','Number_passed','Number_failed','Fraction_passed']) + '\n')
        for eventtype in summarized:
            testpassed = summarized[eventtype]['test'][0]
            testfailed = summarized[eventtype]['test'][1]
            ctrlpassed = summarized[eventtype]['ctrl'][0]
            ctrlfailed = summarized[eventtype]['ctrl'][1]
            testfracpassed = testpassed / float(testpassed + testfailed)
            ctrlfracpassed = ctrlpassed / float(ctrlpassed + ctrlfailed)
            outfh.write(('\t').join([eventtype, 'test', str(testpassed), str(testfailed), str(testfracpassed)]) + '\n')
            outfh.write(('\t').join([eventtype, 'ctrl', str(ctrlpassed), str(ctrlfailed), str(ctrlfracpassed)]) + '\n')
        outfh.close()
