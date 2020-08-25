#This is made for ALEs with 2 isoforms.  The idea is to take each isoform independently.  Using a gtf that contains
#gene names in the gene_name field, intersect the coordinates for each isoform of the event with each entry in the 
#gtf.  Only consider as "hits" gene_names from the gtf file that overlap both isoforms.  If more than one gene_name
#satisfies this, an error is printed and the gene_name is not added to the output.

import pysam
import gffutils
import os
import argparse

def getEvents(eventsfile):
    events = []
    eventsfh = open(eventsfile, 'r')
    for line in eventsfh:
        line = line.strip()
        events.append(line)
    eventsfh.close()
    print 'Entered {0} events.'.format(len(events))
    return events

def eventstocoords(events, annotation):
    gff_fn = annotation
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, force=True)

    db = gffutils.FeatureDB(db_fn)

    coords = {} # {event : [[isoform1chrm, isoform1start, isoform1end, isoform1strand], [isoform2chrm, isoform2start, isoform2end, isoform2strand]]}

    annotatedevents = db.features_of_type('gene')

    for annotatedevent in annotatedevents:
        for event in events:
            eventID = annotatedevent.id
            isoformcoords = []
            if event == eventID:
                for isoform in db.children(annotatedevent, featuretype = 'mRNA'):
                    isoformcoords.append([str(isoform.chrom), isoform.start, isoform.end, str(isoform.strand)])
                coords[event] = isoformcoords

    print 'Found {0} of {1} events in the annotation.'.format(len(coords), len(events))
    for event in events:
        if event not in coords:
            print event
    
    #Cleanup
    os.remove(db_fn)
    if os.path.isfile(db_fn + '-shm'):
        os.remove(db_fn + '-shm')
    if os.path.isfile(db_fn + '-wal'):
        os.remove(db_fn + '-wal')
       
    return coords

def getgenenames(coords, gtf):
    eventstogenes = {} # {event : gene}
    tabixfile = pysam.Tabixfile(gtf)
    for event in coords:
        isoform1genes = []
        isoform2genes = []
        isoform1chrm = coords[event][0][0]
        isoform1start = coords[event][0][1]
        isoform1end = coords[event][0][2]
        isoform1strand = coords[event][0][3]
        isoform2chrm = coords[event][1][0]
        isoform2start = coords[event][1][1]
        isoform2end = coords[event][1][2]
        isoform2strand = coords[event][1][3]
        for entry in tabixfile.fetch(isoform1chrm, isoform1start, isoform1end, isoform1strand, parser = pysam.asGTF()):
            isoform1genes.append(entry.gene_name)
        for entry in tabixfile.fetch(isoform2chrm, isoform2start, isoform2end, isoform2strand, parser = pysam.asGTF()):
            isoform2genes.append(entry.gene_name)

        #Collapse all duplicates
        isoform1genes = list(set(isoform1genes))
        isoform2genes = list(set(isoform2genes))
        #Get genes that overlap both isoforms
        isoformintersection = set(isoform1genes).intersection( set(isoform2genes) )
        
        if len(isoformintersection) > 1:
            print 'WARNING: more than one gene found for event {0}.'.format(event)    
            print event, list(isoformintersection)
        elif len(isoformintersection) == 1:
            eventstogenes[event] = list(isoformintersection)[0]
        elif len(isoformintersection) == 0:
            print 'No gene found for event {0}!!'.format(event)

    print 'Found genes for {0} of {1} events.'.format(len(eventstogenes), len(coords))
    return eventstogenes

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type = str, help = 'List of MISO events you are interested in.', required = True)
    parser.add_argument('--annotations', type = str, help = 'MISO annotations containing events you are interested in.', required = True)
    parser.add_argument('--gtf', type = str, help = 'Gtf annotation of genome. Ones from cuffmerge work well. Must be compressed and tabix indexed with the index in the same directory as the gtf.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    events = getEvents(args.events)
    coords = eventstocoords(events, args.annotations)
    eventstogenes = getgenenames(coords, args.gtf)
    outfh = open(args.outfile, 'w')
    for event in eventstogenes:
        outfh.write(str(event) + '\t' + str(eventstogenes[event]) + '\n')
    outfh.close()
    
