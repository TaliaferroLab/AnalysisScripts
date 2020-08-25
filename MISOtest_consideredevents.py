#Given a directory of misodirectories (each of a different sample), get the events that were
#analyzed in every sample.  These are taken from the batch-logs files.  These events may not have
#passed the initial read count filter, but at least they were looked at.

#Then, given that list of common events, take bed files produced by parseMISO scripts 
#(that are also in the directory of misodirectories) and filter them to only contain events
#that are in the list of common events. 

import os
import sys
import argparse

def getAnalyzedEvents(misodir):
    #Get the events that were looked at in a given miso output directory
    analyzedevents = []
    batchlogs = [os.path.join(os.path.abspath(misodir), 'batch-logs', f) for f in os.listdir(os.path.join(os.path.abspath(misodir), 'batch-logs'))]
    for batchlog in batchlogs:
        infh = open(batchlog, 'r')
        for line in infh:
            line = line.strip()
            if line == 'Computing Psi for 1 genes...':
                eventline = next(infh).strip()
                event = eventline.split('- ')[1]
                analyzedevents.append(event)
        infh.close()
    print '{0} analyzed events in {1}.'.format(len(analyzedevents), os.path.basename(misodir))
    return analyzedevents

def getCommonEvents(dirofmisodirs):
    #Get the events common to a list of lists of analyzed events
    allanalyzedevents = [] #list of lists of analyzed events in each misodir
    for misodir in os.listdir(dirofmisodirs):
        if os.path.isdir(misodir):
            analyzedevents = getAnalyzedEvents(misodir)
            allanalyzedevents.append(analyzedevents)
    commonevents = list(set.intersection(*map(set, allanalyzedevents)))
    print '{0} events in common amongst all samples.'.format(len(commonevents))
    return commonevents

def filterbedfiles(dirofmisodirs, commonevents):
    bedfiles = [f for f in os.listdir(os.path.abspath(dirofmisodirs)) if os.path.isfile(f) and f.endswith('.bed')]
    for bedfile in bedfiles:
        eventcounter = 0
        commoneventcounter = 0
        infh = open(bedfile, 'r')
        outfh = open(bedfile[:-3] + 'filtered.bed', 'w')
        for line in infh:
            eventcounter +=1
            line = line.strip().split('\t')
            eventname = line[3]
            for commonevent in commonevents:
                #if an event in commonevents can be found within eventname
                #The 'in' part is necessary because sometimes parseMISO adds a little bit to the event name
                if commonevent in eventname: 
                    commoneventcounter +=1
                    outfh.write(('\t').join(line) + '\n')
                    break
        infh.close()
        outfh.close()
        print '{0} of {1} events for {2} were events common to all samples.'.format(commoneventcounter,eventcounter,bedfile) 
    
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirofmisodirs', type = str, help = 'Directory containing one MISO output directory per sample.  It should also contain one bed file per MISO directory.  This bed file was produced using parseMISO scripts.')
    args = parser.parse_args()

    numberofmisodirectories = 0
    numberofbedfiles = 0
    
    for f in os.listdir(args.dirofmisodirs):
        if os.path.isdir(f):
            numberofmisodirectories +=1
        elif os.path.isfile(f) and f.endswith('.bed'):
            numberofbedfiles +=1

    if numberofmisodirectories != numberofbedfiles:
        print 'Error! Number of MISO summary bed files in this directory does not equal the number of MISO output directories in this directory!!'
        sys.exit()

    commonevents = getCommonEvents(args.dirofmisodirs)
    filterbedfiles(args.dirofmisodirs, commonevents)

    


