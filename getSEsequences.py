#Takes a list of MISO SE events of interest and the annotation they came from. Returns fasta files of those SE, 200 nt upstream, and 200 nt downstream.

#Usage: python getSEsequences.py --help

import sys
import os
import argparse
import random
from Bio import SeqIO
import gffutils
import gzip

def get_control_events(MISOevents, MISOannotations, numberofctrls):
    #Returns list of events randomly chosen from MISOannotations, provided they are NOT in MISOevents.
    #Returns numberofctrls events.
    events = []
    possiblectrlevents = []
    eventsfh = open(MISOevents, 'r')
    for line in eventsfh:
        line = line.strip()
        events.append(line)
    eventsfh.close()

    #Make gff database of MISOannotations.
    gff_fn = MISOannotations
    db_fn = os.path.basename(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, merge_strategy = 'create_unique')

    db = gffutils.FeatureDB(db_fn)

    annotatedevents = db.features_of_type('gene') #each miso event is of the 'gene' level
    for annotatedevent in annotatedevents:
        if str(annotatedevent.attributes['ID'][0]) not in events:
            possiblectrlevents.append(str(annotatedevent.attributes['ID'][0]))
                                      
    controlevents = random.sample(possiblectrlevents, numberofctrls)
    
    return controlevents
    
def get_regions_of_interest(MISOevents, MISOannotations, controlmode):
    #Returns dictionary of SE coords, upstream 200 nt coords, downstream 200 nt coords
    #MISOevents is a list of events of interest
    exoncoords = []
    if controlmode == True:
        events = MISOevents
    elif controlmode == False:
        events = []
        eventsfh = open(MISOevents, 'r')
        for line in eventsfh:
            line = line.strip()
            events.append(line)
        eventsfh.close()

    regionsofinterest = {}

    #Make gff database of MISOannotations
    gff_fn = MISOannotations
    db_fn = os.path.basename(gff_fn) + '.db'
    
    if os.path.isfile(db_fn) == False:
        print 'Creating annotation database...'
        gffutils.create_db(gff_fn, db_fn, merge_strategy = 'create_unique', verbose = True)

    db = gffutils.FeatureDB(db_fn)

    annotatedevents = db.features_of_type('gene') #each miso event is of the 'gene' level
    for annotatedevent in annotatedevents:
        for event in events:
            if event == str(annotatedevent.attributes['ID'][0]): #if you've found annotated event that matches
                for exon in db.children(annotatedevent, level = 2, featuretype = 'exon'):
                    if exon.attributes['ID'][0].endswith('.se') or exon.attributes['ID'][0].endswith('.sk'): #if this is the skipped exon
                        chrm = annotatedevent.chrom
                        strand = annotatedevent.strand
                        SEcoords = [exon.start, exon.stop]
                        if strand == '+':
                            upstreamcoords = [int(exon.start) - 200, int(exon.start) - 1]
                            downstreamcoords = [int(exon.stop) + 1, int(exon.stop) + 200]
                        elif strand == '-':
                            upstreamcoords = [int(exon.stop) + 1, int(exon.stop) + 200]
                            downstreamcoords = [int(exon.start) - 200, int(exon.start) - 1]

                        regionsofinterest[event] = [{'chrm':chrm}, {'strand':strand}, {'SEcoords': SEcoords}, 
                                                    {'UPcoords':upstreamcoords}, {'DNcoords':downstreamcoords}]

                        exoncoords.append([chrm, 'SE', 'dn_intron', str(downstreamcoords[0]), str(downstreamcoords[1]), '.', strand, '.', 'ID=' + event])

    #outfh = open('ControlIncluded_downstream.gff', 'w')
    #for event in exoncoords:
    #	outfh.write(('\t').join(event) + '\n')
    #outfh.close()

    sys.stderr.write('Found annotations for {0} of {1} events.\n'.format(len(regionsofinterest), len(events)))
    
    os.remove(db_fn)
    #os.remove(db_fn + '-shm')
    #os.remove(db_fn + '-wal')
    return regionsofinterest

def getSequences(genomefasta, regionsofinterest):
    #Takes directory of chrom sequences and regionsofinterest dict from above and returns sequences for each region

    #seqdirectory = os.path.abspath(DirOfChrmSeqs)
    #seqfiles = [os.path.join(seqdirectory, file) for file in os.listdir(seqdirectory)]

    #Make index of fasta files
    sys.stderr.write('Indexing chromosomes...\n')
    seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
    sys.stderr.write('{0} chromosomes indexed.\n'.format(len(seq_dict)))
    
    SEfasta = {}
    Upstreamfasta = {}
    Downstreamfasta = {}
    Combinedfasta = {}

    for event in regionsofinterest:
        for attribute in regionsofinterest[event]:
            if attribute.has_key('chrm'):
                chrm = attribute.values()[0]
            if attribute.has_key('strand'):
                strand = attribute.values()[0]
            if attribute.has_key('SEcoords'):
                SEcoords = attribute.values()[0]
            if attribute.has_key('UPcoords'):
                UPcoords = attribute.values()[0]
            if attribute.has_key('DNcoords'):
                DNcoords = attribute.values()[0]

        fastaID = '>' + event
        if strand == '+':
            SEseq = seq_dict[chrm].seq[SEcoords[0]-1:SEcoords[1]].upper()
            Upstreamseq = seq_dict[chrm].seq[UPcoords[0]-1:UPcoords[1]].upper()
            Downstreamseq = seq_dict[chrm].seq[DNcoords[0]-1:DNcoords[1]].upper()
            Combinedseq = Upstreamseq + SEseq + Downstreamseq
        elif strand == '-':
            SEseq = seq_dict[chrm].seq[SEcoords[0]-1:SEcoords[1]].upper().reverse_complement()
            Upstreamseq = seq_dict[chrm].seq[UPcoords[0]-1:UPcoords[1]].upper().reverse_complement()
            Downstreamseq = seq_dict[chrm].seq[DNcoords[0]-1:DNcoords[1]].upper().reverse_complement()
            Combinedseq = Upstreamseq + SEseq + Downstreamseq
        SEfasta[fastaID] = SEseq
        Upstreamfasta[fastaID] = Upstreamseq
        Downstreamfasta[fastaID] = Downstreamseq
        Combinedfasta[fastaID] = Combinedseq

    
    #os.remove('seqdb.idx')
    
    return SEfasta, Upstreamfasta, Downstreamfasta, Combinedfasta
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--misoevents', type = str, help = 'list of events of interest')
    parser.add_argument('--misoannotations', type = str, help = 'MISO annotation file containing events of interest and possibly others')
    parser.add_argument('--controls', type = int, help = 'Optional. Pick <int> control events from misoannotations that are NOT in misoevents. Retrieve their sequences. If chosen, will return control sequences, NOT sequences from misoevents file.')
    parser.add_argument('--genomefasta', type = str, help = 'Fasta file of genome.')
    parser.add_argument('--SEoutput', type = str, help = 'Output for SE sequence file.')
    parser.add_argument('--UPoutput', type = str, help = 'Output for upstream sequence file.')
    parser.add_argument('--DNoutput', type = str, help = 'Output for downstream sequence file.')
    parser.add_argument('--Combinedoutput', type = str, help = 'Output for combined sequence file.') 
    args = parser.parse_args()

    if args.controls:
        sys.stderr.write('Operating in control mode!! Retrieving sequences for {0} control events.\n'.format(args.controls))
        controlevents = get_control_events(args.misoevents, args.misoannotations, args.controls)
        regionsofinterest = get_regions_of_interest(controlevents, args.misoannotations, True)
    else:
        regionsofinterest = get_regions_of_interest(args.misoevents, args.misoannotations, False)

    SEfasta, Upstreamfasta, Downstreamfasta, Combinedfasta = getSequences(args.genomefasta, regionsofinterest)
    
    SEfh = open(args.SEoutput, 'w')
    UPfh = open(args.UPoutput, 'w')
    DNfh = open(args.DNoutput, 'w')
    Combinedfh = open(args.Combinedoutput, 'w')

    for entry in SEfasta:
        SEfh.write(entry + '\n' + str(SEfasta[entry]) + '\n')
    for entry in Upstreamfasta:
        UPfh.write(entry + '\n' + str(Upstreamfasta[entry]) + '\n')
    for entry in Downstreamfasta:
        DNfh.write(entry + '\n' + str(Downstreamfasta[entry]) + '\n')
    for entry in Combinedfasta:
        Combinedfh.write(entry + '\n' + str(Combinedfasta[entry]) + '\n')

    SEfh.close()
    UPfh.close()
    DNfh.close()
    Combinedfh.close()
