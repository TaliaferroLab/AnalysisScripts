#This script is designed to take a PSI table produced by MISOSigEventsTimecourse.py
#and assign genes to events in that table. It does this in three steps.  
#First it looks to see if the event boundaries lie completely within the boundaries of one
#and only one known transcript. If it doesn't, it tries to get an ENS ID from the gene to
#ENS ID files provided with MISO annotations. If that doesn't work, then it tries to derive
#and ENS ID from the event name.  Sometimes event names contains ENS IDs, especially those
#produced by Jason. If that also fails, then the gene name is defined as 'unknown.'
#It then converts any ENS ID that it retrieved to a gene short name.

#The typical transcript gtf to use is EnsemblTranscripts.mm9.gtf.  This was retrieved from
#UCSC table browser. The conversion of ENS IDs to gene short names is done using 
#Ensembl2ShortGeneName.txt.  The final output is a reproduced PSI table with short gene names woven in.

#Usage: python AddGenesNamestoPSItable.py -h

import os
import sys
import argparse

def getEvents(PSItable):
    #Get all events from a PSItable file and organize them by eventtype
    eventdict = {} # {eventtype : [events]}
    with open(PSItable, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            eventname = line[0]
            eventtype = line[1]
            if 'Event' not in eventname: #skip header
                if eventtype not in eventdict:
                    eventdict[eventtype] = [eventname]
                elif eventtype in eventdict:
                    eventdict[eventtype].append(eventname)

    return eventdict

def getEventCoords(eventtypes, MISOgffs, eventdict):
    #Get 'gene' boundaries for all events of the eventtypes of interest
    #If an event is one of the ones in the PSI table (i.e. it's in eventdict)
    #then record its chrm, start, stop, and strand.
    gfffiledict = dict(zip(eventtypes, MISOgffs)) # {eventtype1 : MISOgffile1, eventtype2 : MISOgffile2}
    coorddict = {} # {eventtype : {eventname : [chrm, start, stop, strand]}}
    counter = 0
    for eventtype in eventtypes:
        coorddict[eventtype] = {}
        gfffile = open(os.path.expanduser(gfffiledict[eventtype]), 'r') #This turns ~ into the home directory
        for line in gfffile:
            counter +=1
            if counter % 10000 == 0:
                print counter
            line = line.strip().split('\t')
            if len(line) > 1 and line[2] == 'gene': #ignore comment lines
                chrm = line[0]
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                eventname = line[8].split(';')[0][3:]
                if eventname in eventdict[eventtype]:
                    #Worried about some events around ends of genes may have boundaries outside known
                    #transcript boundaries. Chopping 100 nt off the most distal coordinates of each event
                    #to try to mitigate this.
                    if eventtype == 'AFE' and strand == '+':
                        coorddict[eventtype][eventname] = [chrm, start + 100, stop, strand]
                    elif eventtype == 'AFE' and strand == '-':
                        coorddict[eventtype][eventname] = [chrm, start, stop - 100, strand]
                    elif eventtype == 'ALE' or eventtype == 'TandemUTR' and strand == '+':
                        coorddict[eventtype][eventname] = [chrm, start, stop - 100, strand]
                    elif eventtype == 'ALE' or eventtype == 'TandemUTR' and strand == '-':
                        coorddict[eventtype][eventname] = [chrm, start + 100, stop, strand]
                    else:
                        coorddict[eventtype][eventname] = [chrm, start, stop, strand]

        gfffile.close()

    for eventtype in coorddict:
        print 'Retrieved {0} {1} events from MISO annotations.'.format(eventtype, len(coorddict[eventtype]))
    
    return coorddict

def parsegtf(gtf):
    transcriptdict = {} #{transcriptID : [chrm, start, stop, strand, geneID]}
    with open(gtf, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0].startswith('ENS'): #skip header
                transcriptID = line[0]
                chrm = line[1]
                strand = line[2]
                start = int(line[3])
                stop = int(line[4])
                geneID = line[7]
                transcriptdict[transcriptID] = [chrm, start, stop, strand, geneID]
                
    print 'Parsed info on {0} transcripts.'.format(len(transcriptdict))
    return transcriptdict

def parsegivengenes(eventtypes, eventstogenesfiles):
    #Parse the MISO-supplied events to ENSgene ID files.
    MISOprovidedeventstogenes = {} #{eventtype : {eventname : ENSgeneID}}
    filedict = dict(zip(eventtypes, eventstogenesfiles)) # {eventtype1 : ENSgenefile1, eventtype2 : ENSgenefile2}
    for eventtype in filedict:
        MISOprovidedeventstogenes[eventtype] = {}
        ENSfile = open(os.path.expanduser(filedict[eventtype]), 'r')
        for line in ENSfile:
            line = line.strip().split('\t')
            if line[0] != 'event_id' and ',' not in line[1] and line[1].startswith('ENSMUS'):
                eventname = line[0]
                gene = line[1]
                MISOprovidedeventstogenes[eventtype][eventname] = gene

        ENSfile.close()

    return MISOprovidedeventstogenes

def intersectcoords(coorddict, transcriptdict, MISOprovidedeventstogenes):
    eventstogenes = {} #{eventtype : {eventname : geneID}}
    for eventtype in coorddict:
        eventstogenes[eventtype] = {}
        unknownevents = 0
        counter = 0
        for eventname in coorddict[eventtype]:
            counter +=1
            if counter % 10000 == 0:
                print 'Intersecting {0} event {1} of {2}.'.format(eventtype, counter, len(coorddict[eventtype]))
            #print 'Intersecting {0} event {1}...'.format(eventtype, eventname)
            geneintersections = []
            eventchrm = coorddict[eventtype][eventname][0]
            eventstart =  coorddict[eventtype][eventname][1]
            eventstop =  coorddict[eventtype][eventname][2]
            eventstrand =  coorddict[eventtype][eventname][3]
            for transcript in transcriptdict:
                transcriptchrm = transcriptdict[transcript][0]
                transcriptstart = transcriptdict[transcript][1]
                transcriptstop = transcriptdict[transcript][2]
                transcriptstrand = transcriptdict[transcript][3]
                geneID = transcriptdict[transcript][4]
                if eventchrm == transcriptchrm and eventstart >= transcriptstart and eventstop <= transcriptstop and eventstrand == transcriptstrand:
                    geneintersections.append(geneID)

            #remove duplicate geneIDs
            geneintersections = list(set(geneintersections))

            if len(geneintersections) == 1:
                eventstogenes[eventtype][eventname] = geneintersections[0]
            elif len(geneintersections) == 0 or len(geneintersections) > 1:
                if eventname in MISOprovidedeventstogenes[eventtype]:
                    eventstogenes[eventtype][eventname] = MISOprovidedeventstogenes[eventtype][eventname]
                elif eventname.startswith('ENSMUSG'):
                    eventstogenes[eventtype][eventname] = eventname.split('@')[0]
                elif eventname.split(':')[-1].startswith('ENSMUSG'):
                    eventstogenes[eventtype][eventname] = eventname.split(':')[-1]
                elif '@' not in eventname and ':' in eventname and eventtype == 'TandemUTR':
                    eventstogenes[eventtype][eventname] = eventname.split(':')[1]
                elif '@' not in eventname:
                    eventstogenes[eventtype][eventname] = eventname
                else:
                    #print 'No gene found for {0} event {1}!'.format(eventtype, eventname)
                    eventstogenes[eventtype][eventname] = 'unknown'
                    unknownevents +=1

            

        print 'Could not find ensembl IDs for {0} of {1} {2} events.'.format(unknownevents, len(coorddict[eventtype]), eventtype)

    return eventstogenes

def ensemblID2shortname(ensembltoshort, eventstogenes):
    #Cross reference the ensembl gene IDs with a tab separated list of ensembl IDs and gene short names (ensembltoshort)
    ensembl2shortdict = {}
    events2IDs = [] # [eventname, eventtype, ensemblID, shortname]
    with open(ensembltoshort, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            ensemblid = line[0]
            shortname = line[1]
            ensembl2shortdict[ensemblid] = shortname

    for eventtype in eventstogenes:
        geneswithIDs = 0
        for eventname in eventstogenes[eventtype]:
            if eventstogenes[eventtype][eventname] in ensembl2shortdict:
                geneID = eventstogenes[eventtype][eventname]
                shortname = ensembl2shortdict[geneID]
                geneswithIDs +=1
                events2IDs.append([eventname, eventtype, geneID, shortname])
            else:
                geneID = eventstogenes[eventtype][eventname]
                events2IDs.append([eventname, eventtype, geneID, geneID])

        print 'Found short name IDs for {0} of {1} {2} events.'.format(geneswithIDs, len(eventstogenes[eventtype]), eventtype)

    return events2IDs

                
                
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PSItable', type = str, required = True, help = 'PSI table produced by MISOSigEventsTimecoursev2.0.py. This is a table of all events and their PSI values in the various samples.')
    parser.add_argument('--eventtypes', type = str, required = True, help = 'Comma separated list of event types you want to consider.  You must have MISO annotation gff files for each.')
    parser.add_argument('--MISOgffs', type = str, required = True, help = 'Comma separated list of MISO annotation gff files for the event types you specified. MUST BE IN THE SAME ORDER AS THE EVENT TYPES.')
    parser.add_argument('--eventstogenes', type = str, required = True, help = 'Comma separated list of supplied MISO events to genes files in same order as event types.')
    parser.add_argument('--transcriptgtf', type = str, required = True, help = 'Gtf of transcript coords. EnsemblTranscripts.mm9.gtf as an example.')
    parser.add_argument('--ensembltoshort', type = str, required = True, help = 'File of tab separated values of Ensembl Gene ID and associated gene short name.')
    parser.add_argument('--output', type = str, required = True, help = 'Output file.')
    args = parser.parse_args()

    eventtypes = args.eventtypes.split(',')
    MISOgffs = args.MISOgffs.split(',')
    eventstogenesfiles = args.eventstogenes.split(',')

    if len(eventtypes) != len(MISOgffs):
        print 'ERROR: You must specify one MISO gff for each event type.'
        sys.exit()

    if len(eventtypes) != len(eventstogenesfiles):
        print 'ERROR: You must specify one eventstogenes file for each event type.'
        sys.exit()

    eventdict = getEvents(args.PSItable)
    coorddict = getEventCoords(eventtypes, MISOgffs, eventdict)
    transcriptdict = parsegtf(args.transcriptgtf)
    MISOprovidedeventstogenes = parsegivengenes(eventtypes, eventstogenesfiles)
    eventstogenes = intersectcoords(coorddict, transcriptdict, MISOprovidedeventstogenes)
    events2IDs = ensemblID2shortname(args.ensembltoshort, eventstogenes)
    outfh = open(args.output, 'w')
    infh = open(args.PSItable, 'r')
    for line in infh:
        line = line.strip().split('\t')
        if line[0] == 'Event': #if header
            line.insert(1, 'Gene')
            outfh.write(('\t').join(line) + '\n')
        elif line[0] != 'Event':
            eventname = line[0]
            eventtype = line[1]
            for event in events2IDs:
                if eventname == event[0] and eventtype == event[1]:
                    shortname = event[3]
                    line.insert(1, shortname)
                    outfh.write(('\t').join(line) + '\n')

    outfh.close()
    infh.close()
    
