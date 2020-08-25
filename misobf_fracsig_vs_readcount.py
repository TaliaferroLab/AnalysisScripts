#Want to know what fraction of events from a MISO comparison are significant (to some filter).
#Then, bin events by number of informative reads ((0,1 or 1,0 reads)).
#
#Actually, 

import sys
import operator
import numpy as np

def geteventdict(miso_bffile):
    bffh = open(miso_bffile, 'r')
    eventdict = {} # {eventname : [number of informative reads ((1,0) or (0,1)), filter pass (yes/no)]}
    eventcounter = 0
    for line in bffh:
        line = line.strip().split('\t')
        if line[0] == 'event_name': #skip header
            continue
        numberofisoforms = len(line[1].split(',')) + 1
        if numberofisoforms == 2:
            eventcounter +=1
            eventname = line[0]
            filterpass = False
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
            informativereadcount = combinedsample1counts + combinedsample2counts
            
            if bayesfactor >= 10:
                filterpass = True
            eventdict[eventname] = [informativereadcount, deltapsi, filterpass]

    print 'Considered {0} two isoform events.'.format(eventcounter)
    bffh.close()

    return eventdict

def summarizeeventdict(eventdict1, eventdict2):
    #eventdict = {eventname : [number of informative reads ((1,0) or (0,1)), filter pass (yes/no)]}
    commoneventdict = {} # {commonevent : [avg number of informative reads, filter pass in both samples]}

    for event in eventdict1:
        if event in eventdict2:
            commoneventdict[event] = []

    for event in commoneventdict:
        sample1readcount = eventdict1[event][0]
        sample2readcount = eventdict2[event][0]
        avgreadcount = (sample1readcount + sample2readcount) / 2.0
        filterpassinboth = False
        sample1deltapsi = eventdict1[event][1]
        sample2deltapsi = eventdict2[event][1]
        sample1filterpass = eventdict1[event][2]
        sample2filterpass = eventdict2[event][2]
        #Check and see if it passed filters in both samples and also that their delta psi values have the same sign
        if sample1filterpass == True and sample2filterpass == True and sample1deltapsi * sample2deltapsi > 0:
            filterpassinboth = True
        commoneventdict[event] = [avgreadcount, filterpassinboth]

    print '{0} two isoform events present in both samples.'.format(len(commoneventdict))
    
    bindict = {} # {binnumber : {eventname : {readcount : filterpass}}}
    outdict = {} # {binnumber : [avg. read count in bin, percent events that filterpass]}
    numberofevents = len(commoneventdict)
    numberofbins = 20
    eventsperbin = numberofevents / numberofbins
    
    for binnumber in range(1, numberofbins + 1):
        bindict[binnumber] = {}
        lowerbinbound = (binnumber - 1) * eventsperbin
        upperbinbound = (binnumber) * eventsperbin 
        for idx, event in enumerate(sorted(commoneventdict.items(), key = operator.itemgetter(1))):
            #print idx, event, binnumber, lowerbinbound, upperbinbound
            eventname = event[0]
            readcount = event[1][0]
            filterpass = event[1][1]
            if idx + 1 > lowerbinbound and idx + 1 <= upperbinbound:
                bindict[binnumber][eventname] = {}
                bindict[binnumber][eventname][readcount] = filterpass

    
    for binnumber in bindict:
        outdict[binnumber] = []
        readcounts = []
        filterpass = 0
        filterfail = 0
        for event in bindict[binnumber]:
            readcount = bindict[binnumber][event].keys()[0]
            filterresult = bindict[binnumber][event].values()[0]
            readcounts.append(readcount)
            if filterresult == True:
                filterpass +=1
            elif filterresult == False:
                filterfail +=1

        avgreadcount = np.mean(readcounts)
        fractionpass = float(filterpass / float(filterpass + filterfail))
        outdict[binnumber] = [avgreadcount, fractionpass]

    return outdict
    

eventdict1 = geteventdict(sys.argv[1])
eventdict2 = geteventdict(sys.argv[2])
summarizeddict = summarizeeventdict(eventdict1, eventdict2)
outfh = open(sys.argv[4], 'w')
for binnumber in summarizeddict:
    outfh.write(str(binnumber) + '\t' + str(summarizeddict[binnumber][0]) + '\t' + str(summarizeddict[binnumber][1]) + '\t' + sys.argv[3] + '\n')

outfh.close()
