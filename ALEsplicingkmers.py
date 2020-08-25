#Not yet finished.

#Takes a list of ALE events and counts kmers 200 nt upstream and 200 nt downstream of both the proximal and 
#distal ALE 3' splice sites

import os
import gffutils
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np

def getRegions(ALEevents, ALEannotation):
    eventsofinterest = []
    foundevents = 0
    eventsofinterestfh = open(ALEevents, 'r')
    for line in eventsofinterestfh:
        line = line.strip()
        eventsofinterest.append(line)
    eventsofinterestfh.close()
    proximalisoforms = []  #200 nt around proximal 3' splice site of ALE, gff format
    distalisoforms = [] #200 nt around distal 3' splice site of ALE, gff format

    gff_fn = ALEannotation
    db_fn = os.path.basename(gff_fn) + '.db'

    print 'Indexing annotation...'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, force = False, verbose = False)
    db = gffutils.FeatureDB(db_fn)
    events = db.features_of_type('gene')

    for event in events:
        if event.id in eventsofinterest:
            foundevents +=1 #this event exists in the annotation
            if event.strand == '+':
                isoformcounter = 1 
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'end'): #order by end of isoform, this puts proximal isoform first
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter == 1:
                        proximalisoforms.append([isoform.chrom, 'ALEprox3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])
                        isoformcounter +=1
                    elif isoformcounter == 2:
                        distalisoforms.append([isoform.chrom, 'ALEdist3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])

            elif event.strand == '-':
                isoformcounter = 1 
                #order by end of isoform, this puts proximal isoform first
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'start', reverse = True): 
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter == 1:
                        proximalisoforms.append([isoform.chrom, 'ALEprox3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])
                        isoformcounter +=1
                    elif isoformcounter == 2:
                        distalisoforms.append([isoform.chrom, 'ALEdist3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])

    print 'Found {0} of {1} provided events in the annotation.'.format(foundevents, len(eventsofinterest))
    #Cleanup
    os.remove(db_fn)
    return proximalisoforms, distalisoforms

def gfftofasta(gff, genomesequence):
    #Given gff in list form (output of above function, each line is a list), retrieve sequences and return in fasta format
    #gff is list of lists, each line is a single list

    GCs = []
    seqs = {} #dictionary where key is ID and value is sequence
    
    print 'Indexing genome sequence...'
    seq_dict = SeqIO.to_dict(SeqIO.parse(genomesequence, 'fasta'))
    print 'Retrieving sequences...'
    
    for entry in gff:
        chrm = entry[0]
        start = int(entry[3])
        stop = int(entry[4])
        strand = entry[6]
        ID = str(entry[8])

        if strand == '+':
            sequence = str(seq_dict[chrm].seq[start-200:start+200].upper())
        elif strand == '-':
            sequence = str(seq_dict[chrm].seq[stop-200:stop+200].reverse_complement().upper())
        
        seqs[ID] = sequence
        GCs.append(float(GC(sequence)))

    GCcontent = np.mean(GCs)

    return seqs, GCcontent

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type = str, help = 'List of ALE events to consider.', required = True)
    parser.add_argument('--annotations', type = str, help = 'ALE annotations containing events in MISO gff format.', required = True)
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.', required = True)
    args = parser.parse_args()

    proximalisoforms, distalisoforms = getRegions(args.events, args.annotations)
    proximalseqs, proximalGC = gfftofasta(proximalisoforms, args.genomefasta)
    print 'The average GC content of the proximal sequences is {0}.'.format(proximalGC)
    distalseqs, distalGC = gfftofasta(distalisoforms, args.genomefasta)
    print 'The average GC content of the distal sequences is {0}.'.format(distalGC)

    #Write the sequences to temporary fasta files
    outfh = open('prox.temp.fasta', 'w')
    for entry in proximalseqs:
        outfh.write('>' + entry + '\n' + proximalseqs[entry] + '\n')
    outfh.close()
    outfh = open('dist.temp.fasta', 'w')
    for entry in distalseqs:
        outfh.write('>' + entry + '\n' + distalseqs[entry] + '\n')
    outfh.close()
