#Necessary modules: scipy, numpy, biopython, gffutils, R v2.15, rpy2
#Must also have kmerenrichment.py in your pythonpath

#Take a list of tandem UTRs and calculate kmer enrichments between the distal and proximal
#parts of the tandem UTR.

#python TandemUTRkmerenrichment.py --help

import os
import gffutils
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
import numpy as np
import kmerenrichment
import subsampleGC

def getUTRregions(TUTRevents, TUTRannotation):
    #Given a tandem UTR miso annotation and a list of events you are interested in, retrieve the proximal and distal
    #UTRs separately. The proximal UTR will NOT be contained within the distal UTR.  Rather the distal UTR will begin
    #right after the first polyA site.
    eventsofinterest = []
    proximalisoforms = [] #from stop codon to first polyA
    distalisoforms = [] #from first polyA to second polyA
    foundevents = 0
    eventsofinterestfh = open(TUTRevents, 'r')
    for line in eventsofinterestfh:
        line = line.strip()
        eventsofinterest.append(line)
    eventsofinterestfh.close()
    
    gff_fn = TUTRannotation
    db_fn = os.path.basename(gff_fn) + '.db'

    print 'Indexing annotation...'
    gffutils.create_db(gff_fn, db_fn, force = True, verbose = False)
    db = gffutils.FeatureDB(db_fn)
    events = db.features_of_type('gene')
    
    for event in events:
        if event.id in eventsofinterest:
            foundevents +=1 #this event exists in the annotation
            if event.strand == '+':
                isoformcounter = 1
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'end'): #order by end with shortest one first
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter ==1:
                        proximalisoforms.append([isoform.chrom, 'TandemUTR', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])
                        distalstart = int(isoform.end) + 1 #distal isoform starts right after proximal one ends
                        isoformcounter +=1
                    elif isoformcounter ==2:
                        distalisoforms.append([isoform.chrom, 'TandemUTR', 'mRNA', str(distalstart), str(isoform.end), '.', isoform.strand,
                                               '.', isoform.id + ';Parent=' + parentid])
            
            elif event.strand == '-':
                isoformcounter = 1
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'start', reverse = True): #order by start with shortest one first
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter ==1:
                        proximalisoforms.append([isoform.chrom, 'TandemUTR', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid])
                        distalend = int(isoform.start) - 1 #distal isoform ends right before proximal one starts
                        isoformcounter +=1
                    elif isoformcounter == 2:
                        distalisoforms.append([isoform.chrom, 'TandemUTR', 'mRNA', str(isoform.start), str(distalend), '.', isoform.strand, 
                                               '.', isoform.id + ';Parent=' + parentid])
                    
    print 'Found {0} of {1} provided events in the annotation.'.format(len(eventsofinterest), foundevents)
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
            sequence = str(seq_dict[chrm].seq[start:stop+1].upper())
        elif strand == '-':
            sequence = str(seq_dict[chrm].seq[start-1:stop].reverse_complement().upper())
        #Remove last 50 nt to get rid of polyA signals
        if len(sequence) > 50:
            seqs[ID] = sequence[:-50]
            GCs.append(float(GC(sequence)))

    GCcontent = np.mean(GCs)

    return seqs, GCcontent



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type = str, help = 'List of tandem UTR events to consider.', required = True)
    parser.add_argument('--annotations', type = str, help = 'TandemUTR annotations containing events in MISO gff format.', required = True)
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.', required = True)
    parser.add_argument('-k', type = int, help = 'Length of kmers to look for.', required = True)
    parser.add_argument('--output', type = str, help = 'Output file.', required = True)
    parser.add_argument('--GCmatch', dest = 'GCmatch', help = 'Subsample distal sequences to match GC content for proximal sequences. Default is no subsampling.',
                        required = False, action = 'store_true')
    parser.set_defaults(GCmatch = False)
    args = parser.parse_args()

    proximalisoforms, distalisoforms = getUTRregions(args.events, args.annotations)
    outfh = open('distalUTRsections.gff3', 'w')
    for line in distalisoforms:
        outfh.write(('\t').join(line) + '\n')
    outfh.close()
    proximalseqs, proximalGC = gfftofasta(proximalisoforms, args.genomefasta)
    print 'The average GC content of the proximal UTRs is {0}.'.format(proximalGC)
    distalseqs, distalGC = gfftofasta(distalisoforms, args.genomefasta)
    print 'The average GC content of the distal UTRs is {0}.'.format(distalGC)

    #Write the sequences to temporary fasta files
    outfh = open('prox.temp.fasta', 'w')
    for entry in proximalseqs:
        outfh.write('>' + entry + '\n' + proximalseqs[entry] + '\n')
    outfh.close()
    outfh = open('dist.temp.fasta', 'w')
    for entry in distalseqs:
        outfh.write('>' + entry + '\n' + distalseqs[entry] + '\n')
    outfh.close()

    #Subsample to match GC contents...make distal GCs look like proximal GCs
    if args.GCmatch == True:
        #This returns a list of fasta records. Each record is itself a list where the first item is the ID and the second is the sequence.
        subsampleddistal = subsampleGC.subsampleGC('prox.temp.fasta', 'dist.temp.fasta')
        outfh = open('dist.subsampled.temp.fasta', 'w')
        for entry in subsampleddistal:
            outfh.write(str(entry[0]) + '\n' + str(entry[1]) + '\n')
        outfh.close()

    #Calculate kmer enrichments
    print 'Calculating kmer enrichments...'
    if args.GCmatch == True:
        outputlist = kmerenrichment.countKmers('dist.subsampled.temp.fasta', 'prox.temp.fasta', args.k)
    elif args.GCmatch == False:
        outputlist = kmerenrichment.countKmers('dist.temp.fasta', 'prox.temp.fasta', args.k)
    outfh = open(args.output, 'w')

    outfh.write('kmer' + '\t' + 'distal_count' + '\t' + 'proximal_count' + '\t' + 
                'enrichment' + '\t' + 'fisher\'s exact p' + '\t' + 'BH corrected p''\n')
    for kmer in outputlist:
        outfh.write(('\t').join([str(item) for item in kmer]) + '\n')
    outfh.close()

'''
    #Cleanup temporary files
    os.remove('dist.temp.fasta')
    os.remove('prox.temp.fasta')
    if args.GCmatch == True:
        os.remove('dist.subsampled.temp.fasta')
'''  
