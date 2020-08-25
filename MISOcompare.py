#Takes a particular file structure (described below) and counts number of events that pass particular
#BF, deltaPSI, and read count fiters for each comparison in the file structure.  
#These filters can be easily changed in the code.

#Usage: python MISOcompare.py --help

import sys
import os
import argparse

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
            sample2 = comparison.split('.')[1].split('_')[2]
            comparisontype = ''
            
            if 'Soma' in sample1 and 'Axon' in sample2:
                comparisontype = 'compartment'
            if 'Axon' in sample1 and 'Soma' in sample2:
                comparisontype = 'compartment'
            if 'CADSoma' in sample1 and 'CADSoma' in sample2:
                comparisontype = 'replicate'
            if 'CADAxon' in sample1 and 'CADAxon' in sample2:
                comparisontype = 'replicate'
            if 'N2ASoma' in sample1 and 'N2ASoma' in sample2:
                comparisontype = 'replicate'
            if 'N2AAxon' in sample1 and 'N2AAxon' in sample2:
                comparisontype = 'replicate'
            if 'CADSoma' in sample1 and 'N2ASoma' in sample2:
                comparisontype = 'celltype'
            if 'CADAxon' in sample1 and 'N2AAxon' in sample2:
                comparisontype = 'celltype'
            if 'N2ASoma' in sample1 and 'CADSoma' in sample2:
                comparisontype = 'celltype'
            if 'N2AAxon' in sample1 and 'CADAxon' in sample2:
                comparisontype = 'celltype'
            

            
            #if comparisontype == 'replicate' or comparisontype == 'celltype':
            #comparisontype = 'other'
            
            #Get bayes factor file
            bffile = os.path.join(os.path.abspath(eventtypedir), comparison, 'bayes-factors', comparison + '.miso_bf')
            bffh = open(bffile, 'r')
            for line in bffh:
                line = line.strip().split('\t')
                if line[0] == 'event_name': #skip header
                    continue
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
                    if abs(deltapsi) >= 0.1 and bayesfactor >= 10 and sample1countspass and sample2countspass:
                        sigeventcounter +=1

            if sample1 != 'N2AAxon1' and sample2 != 'N2AAxon1' and sample1 != sample2: #had always removed N2AAxon1
                comparison_results = [sample1, sample2, eventtype, comparisontype, str(eventcounter), str(sigeventcounter/float(eventcounter))]
                outlists.append(comparison_results)
                    
            bffh.close()
            
    return outlists

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--MISOdirectory', type = str, help = 'Main MISO Output directory. Contains subdirectories for all event types.  These subdirectories contains directories called "comparisons" which themselves contain directories for each pairwise comparison (e.g. CADAxon1.psi_vs_CADAxon2.psi) which themselves contain a directory called bayes-factors. Required.')
    parser.add_argument('--outfile', type = str, help = 'Output file. Required.')
    args = parser.parse_args()

    outlists = MISOcompare(args.MISOdirectory)
    outfh = open(args.outfile, 'w')
    outfh.write('Sample1' + '\t' + 'Sample2' + '\t' + 'EventType' + '\t' + 'ComparisonType' + '\t' + 'Events' + '\t' + 'SigEvents' + '\n')
    for entry in outlists:
        outfh.write(('\t').join(entry) + '\n')
    outfh.close()
