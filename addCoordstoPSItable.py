#Created by MT on 12/12/14

#Takes a PSI table produced by MISOSigEventsTimecourse_v2.0.py and MISO annotations and adds coordinates
#from the event to the PSI table.  For AFEs and ALEs, it adds coordinates of the most proximal exon in 
#the event.  From MXEs, it adds coords from the first exon.  From RI, it adds coords of the first exon.
#For SE, it adds coords of the skipped exon, and for tandem UTR it adds coords of the core UTR region.

#Usage: python addCoordstoPSItable.py --help

import argparse
import gffutils
import os

def addCoords(AFEannot, ALEannot, MXEannot, RIannot, SEannot, TUTRannot, PSItable, outfile):
    AFEdb_fn = os.path.basename(AFEannot) + '.db'
    ALEdb_fn = os.path.basename(ALEannot) + '.db'
    MXEdb_fn = os.path.basename(MXEannot) + '.db'
    RIdb_fn = os.path.basename(RIannot) + '.db'
    SEdb_fn = os.path.basename(SEannot) + '.db'
    TUTRdb_fn = os.path.basename(TUTRannot) + '.db'

    print 'Building AFE annotation database...'
    gffutils.create_db(AFEannot, AFEdb_fn)
    print 'Building ALE annotation database...'
    gffutils.create_db(ALEannot, ALEdb_fn)
    print 'Building MXE annotation database...'
    gffutils.create_db(MXEannot, MXEdb_fn)
    print 'Building RI annotation database...'
    gffutils.create_db(RIannot, RIdb_fn)
    print 'Building SE annotation database...'
    gffutils.create_db(SEannot, SEdb_fn)
    print 'Building TandemUTR annotation database...'
    gffutils.create_db(TUTRannot, TUTRdb_fn)

    AFEdb = gffutils.FeatureDB(AFEdb_fn)
    ALEdb = gffutils.FeatureDB(ALEdb_fn)
    MXEdb = gffutils.FeatureDB(MXEdb_fn)
    RIdb = gffutils.FeatureDB(RIdb_fn)
    SEdb = gffutils.FeatureDB(SEdb_fn)
    TUTRdb = gffutils.FeatureDB(TUTRdb_fn)

    print 'Finding coordinates...'

    outfh = open(outfile, 'w')
    PSItablefh = open(PSItable, 'r')
    for line in PSItablefh:
        PSIvalues = []
        line = line.strip().split('\t')
        #get header
        if line[0] == 'Event':
            header = [line[0], line[1], 'chrom', 'start', 'stop']
            for idx, field in enumerate(line):
                if idx >=2:
                    header.append(field)
            outfh.write(('\t').join(header) + '\n')
            number_of_fields = len(line)
            #Number of fields - Event - EventType - SigComparisons
            number_of_samples = len(line) - 3
        #skip header
        if line[0] != 'Event':
            event = line[0]
            eventtype = line[1]
            for idx, field in enumerate(line):
                if idx > 1 and idx < (number_of_fields - 1):
                    PSIvalues.append(field)
                    
            
            if eventtype == 'AFE':
                eventannots = AFEdb.features_of_type('gene')
                for eventannot in eventannots:
                    if event == eventannot.id:
                        if eventannot.strand == '+':
                            for idx, exon in enumerate(AFEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if (idx + 1) == sum(1 for _ in (AFEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most downstream exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                        elif eventannot.strand == '-':
                            for idx, exon in enumerate(AFEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if (idx + 1) == sum(1 for _ in (AFEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most upstream exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')
            

            if eventtype == 'ALE':
                eventannots = ALEdb.features_of_type('gene')
                for eventannot in eventannots:
                    if event == eventannot.id:
                        if eventannot.strand == '+':
                            for idx, exon in enumerate(ALEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if (idx + 1) == sum(1 for _ in (ALEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most upstream exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                        elif eventannot.strand == '-':
                            for idx, exon in enumerate(ALEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if (idx + 1) == sum(1 for _ in (ALEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most downstream exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

            

            if eventtype == 'MXE':
                 eventannots = MXEdb.features_of_type('gene')
                 for eventannot in eventannots:
                     if event == eventannot.id:
                         if eventannot.strand == '+':
                             for idx, exon in enumerate(MXEdb.children(eventannot, featuretype='exon', order_by='start')):
                                 if (idx + 1) == sum(1 for _ in (MXEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most downstream exon in the event
                                     start = str(exon.start)
                                     stop = str(exon.stop)
                                     chrm = exon.chrom
                                     outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                         elif eventannot.strand == '-':
                             for idx, exon in enumerate(MXEdb.children(eventannot, featuretype='exon', order_by='start')):
                                 if (idx + 1) == sum(1 for _ in (MXEdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most upstream exon in the event
                                     start = str(exon.start)
                                     stop = str(exon.stop)
                                     chrm = exon.chrom
                                     outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

            
            if eventtype == 'RI':
                 eventannots = RIdb.features_of_type('gene')
                 for eventannot in eventannots:
                     if event == eventannot.id:
                         if eventannot.strand == '+':
                             for idx, exon in enumerate(RIdb.children(eventannot, featuretype='exon', order_by='start')):
                                 if (idx + 1) == sum(1 for _ in (RIdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most downstream exon in the event
                                     start = str(exon.start)
                                     stop = str(exon.stop)
                                     chrm = exon.chrom
                                     outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                         elif eventannot.strand == '-':
                             for idx, exon in enumerate(RIdb.children(eventannot, featuretype='exon', order_by='start')):
                                 if (idx + 1) == sum(1 for _ in (RIdb.children(eventannot, featuretype='exon', order_by='start'))): #if this is the most upstream exon in the event
                                     start = str(exon.start)
                                     stop = str(exon.stop)
                                     chrm = exon.chrom
                                     outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

            
            if eventtype == 'SE':
                eventannots = SEdb.features_of_type('gene')
                for eventannot in eventannots:
                    if event == eventannot.id:
                        if eventannot.strand == '+':
                            for idx, exon in enumerate(SEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if idx == 2: #if this is the middle (3 of 5) exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                        elif eventannot.strand == '-':
                            for idx, exon in enumerate(SEdb.children(eventannot, featuretype='exon', order_by='start')):
                                if idx == 2: #if this is the middle (3 of 5) exon in the event
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

            
            if eventtype == 'TandemUTR':
                eventannots = TUTRdb.features_of_type('gene')
                for eventannot in eventannots:
                    if event == eventannot.id:
                        if eventannot.strand == '+':
                            for idx, exon in enumerate(TUTRdb.children(eventannot, featuretype='exon', order_by='end')):
                                if idx == 0: #if this exon has the most upstream stop
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')

                        elif eventannot.strand == '-':
                            for idx, exon in enumerate(TUTRdb.children(eventannot, featuretype='exon', order_by='start', reverse=True)):
                                if idx == 0: #if this is the most upstream start
                                    start = str(exon.start)
                                    stop = str(exon.stop)
                                    chrm = exon.chrom
                                    outfh.write(event + '\t' + eventtype + '\t' + chrm + '\t' + start + '\t' + stop + '\t' + ('\t').join(PSIvalues) + '\t' + line[-1] + '\n')


                
    outfh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PSItable', type = str, help = 'PSI table produced by MISOSigEventTimecourse_v2.0.py.')
    parser.add_argument('--AFEannot', type = str, help = 'AFE annotations in MISO format.')
    parser.add_argument('--ALEannot', type = str, help = 'ALE annotations in MISO format.')
    parser.add_argument('--MXEannot', type = str, help = 'MXE annotations in MISO format.')
    parser.add_argument('--RIannot', type = str, help = 'RI annotations in MISO format.')
    parser.add_argument('--SEannot', type = str, help = 'SE annotations in MISO format.')
    parser.add_argument('--TUTRannot', type = str, help = 'TandemUTR annotations in MISO format.')
    
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    addCoords(args.AFEannot, args.ALEannot, args.MXEannot, args.RIannot, args.SEannot, args.TUTRannot, args.PSItable, args.outfile)
