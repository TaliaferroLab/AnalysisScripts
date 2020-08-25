import gffutils
import argparse
import os

def getstopcodons(gff):
    stopcodons = {} # {stopcodon_n : [chrm, start, stop, strand]}
    #GFF must have stop codons explicitly annotated as such in the 3rd field
    counter = 1
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        if line[2] == 'stop_codon':
            stopcodons['stop_codon_{0}'.format(counter)] = [line[0], line[3], line[4], line[6]]
            counter +=1
    print 'Found {0} stop codons in gff.'.format(len(stopcodons))
    gfffh.close()

    return stopcodons

def stopsinSEs(MISOannotations, stopcodons, outfile):
    PCEevents = 0
    currentevent = ''
    gff_fn = MISOannotations
    db_fn = os.path.basename(gff_fn) + '.db'
    print 'Indexing annotation...'
    if os.path.isfile(db_fn) == False:
        #Use keep_order to keep mRNAs (isoforms) in proper order
        gffutils.create_db(gff_fn, db_fn, force = True, keep_order = True, verbose = False)
    print 'Done indexing!'
    db = gffutils.FeatureDB(db_fn)
    events = db.features_of_type('gene')

    for event in events:
        for exon in db.children(event, featuretype = 'exon'):
            if exon.id[-2:] == 'se' or exon.id[-2:] == 'sk': #if this is the skipped exon
                skippedexon = exon
                for stopcodon in stopcodons:
                    stopcodonchrm = stopcodons[stopcodon][0]
                    stopcodonstart = int(stopcodons[stopcodon][1])
                    stopcodonstop = int(stopcodons[stopcodon][2])
                    stopcodonstrand = stopcodons[stopcodon][3]
                    if exon.chrom == stopcodonchrm and exon.start <= stopcodonstart and exon.end >= stopcodonstop \
                        and exon.strand == stopcodonstrand: #if there's a stop codon in the skipped exon
                        PCEevents +=1
                        with open(outfile, 'a') as f:
                            f.write('stopcdn:' + str(stopcodonchrm) + str(':') + str(stopcodonstart) + str('-') + str(stopcodonstop) + stopcodonstrand + '\n')
                            #if we are still on the same event, don't write the event line again
                            if event.id != currentevent:
                                f.write(str(event) + '\n')

                            #if this is a two isoform event
                            if len(list(db.children(event, featuretype = 'mRNA'))) == 2: 
                                for isoform in db.children(event, featuretype = 'mRNA'):
                                    f.write(str(isoform) + '\n')
                                for eventexon in db.children(event, featuretype = 'exon'):
                                    f.write(str(eventexon) + '\n')

                            #if this event has more than two isoforms
                            elif len(list(db.children(event, featuretype = 'mRNA'))) > 2:
                                inclusionisoform = db.parents(skippedexon, featuretype = 'mRNA')
                                for isoform in db.children(event, featuretype = 'mRNA'):
                                    #If this is the exclusion isoform (sister to the inclusion isoform)
                                    if isoform.id[:-2] == inclusionisoform.id:
                                        exclusionisoform = isoform
                                f.write(str(inclusionisoform) + '\n')
                                for inclusionisoformexon in db.children(inclusionisoform, featuretype = 'exon'):
                                    f.write(str(inclusionisoformexon) + '\n')
                                f.write(str(exclusionisoform) + '\n')
                                for exclusionisoformexon in db.children(exclusionisoform, featuretype = 'exon'):
                                    f.write(str(exclusionisoformexon) + '\n')

                            currentevent = event.id
                                
                        #Only care if it's hit at least once, so break after first instance of stop codon in SE
                        break

    print 'Found {0} events with stop codons in the skipped exon.'.format(PCEevents)

    #Cleanup
    if os.path.isfile(db_fn):
        os.remove(db_fn)
    if os.path.isfile(db_fn + '-shm'):
        os.remove(db_fn + '-shm')
    if os.path.isfile(db_fn + 'wal'):
        os.remove(db_fn + 'wal')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stopcodons', type = str, help = 'Gff file containing stop codons annotated as such in 3rd field.  mm9.genes.chr.jess.gff works well.')
    parser.add_argument('--misoannotations', type = str, help = 'MISO annotations of skipped exon events.')
    parser.add_argument('--output', type = str, help = 'Output gff in MISO annotation format.')
    args = parser.parse_args()

    stopcodons = getstopcodons(args.stopcodons)
    stopsinSEs(args.misoannotations, stopcodons, args.output)
    
