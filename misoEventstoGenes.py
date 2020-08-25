#Takes a list of miso events and returns the genes they are in.  Takes as input a tab separated file with events in the first column and event type in the second.  Also takes a tab separated file of event name and gene name to cross-reference to.  This file was included in the MISO annotations.

#Usage: python misoEventstoGenes.py <events.txt> <crossreference.txt> <output.txt>


import sys

conversionfh = open(sys.argv[2], 'r')
eventsfh = open(sys.argv[1],'r')
outfh = open(sys.argv[3], 'w')

events = []
conversiondict = {}
genenames = []


for line in eventsfh:
    line = line.strip()
    line = line.split('\t')
    events.append(line)

for line in conversionfh:
    line = line.strip()
    line = line.split('\t')
    if line[0] != "event_id":
        conversiondict[line[0]] = line[1]

for event in events:
    eventname = event[0]
    eventtype = event[1]
    genename = []

    if eventtype == 'AFE' or eventtype == 'ALE' or eventtype == 'MXE' or eventtype == 'RI':
        if 'ENSMUS' in eventname:
            genename = eventname[0:18]
        elif '@' not in eventname:
            genename = eventname
        elif eventname in conversiondict:
            genename = conversiondict[eventname]
        else:
            print 'WARNING: no gene name found for event ' + eventname + '\n'

    elif eventtype == 'SE':
        if 'ENSMUS' in eventname:
            genename = eventname[-18:]
        elif '@' not in eventname:
            genename = eventname
        elif eventname in conversiondict:
            genename = conversiondict[eventname]
        else:
            print 'WARNING: no gene name found for event ' + eventname + '\n'

    elif eventtype == 'TandemUTR':
        if '@' not in eventname:
            eventname = eventname.split(':')
            genename = eventname[0]
        elif eventname in conversiondict:
            genename = conversiondict[eventname]
        else:
            print 'WARNING: no gene name found for event ' + eventname + '\n'

    genename = [genename, eventtype]

    #If genename[0] is not empty
    if genename[0]:
        genenames.append(genename)



for gene in genenames:
    '''
    outfh.write(str(gene[0]) + '\t' + str(gene[1]) + '\n')
    '''
    outfh.write(str(gene[0] + '\n'))
eventsfh.close()
conversionfh.close()
outfh.close()
