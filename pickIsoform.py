#Picks desired isoforms from a sorted gff.  Presumes gff is sorted with distal isoform first (actually redoes the sorting itself since gffutils puts the gff back in its original order).  Will skip any isoforms for which an isoform with the same [chromosome, start, stop, strand, and combined length of child exons] has already been reported.

#Usage: python pickIsoform.py <gff_fn> <AFE/ALE> <distal/proximal> > output.gff3

import sys
import os
import gffutils
from operator import itemgetter

gff_fn = sys.argv[1]
event_type = sys.argv[2] #'AFE' or 'ALE'
isoform_to_pick = sys.argv[3] #'proximal' or 'distal'
db_fn = os.path.basename(gff_fn) + '.db'

if os.path.isfile(db_fn) == False: #if database doesn't exist, create it
    gffutils.create_db(gff_fn, db_fn)

db = gffutils.FeatureDB(db_fn)

genes = db.features_of_type('gene')
pickedisoforms = []

for gene in genes:
    print gene
    isoformcounter = 0
    bounds = []
    for isoform in db.children(gene, featuretype = 'mRNA'):
        bounds.append([isoform.start, isoform.stop])

    for isoform in db.children(gene, featuretype = 'mRNA'):
        if event_type == 'AFE':
            if gene.strand == '+':
                if bounds[0][0] != bounds[1][0]:
                    bounds = sorted(bounds, key = itemgetter(0)) #sort isoform start/stop by start coord, lowest to highest
                elif bounds[0][0] == bounds[1][0]:
                    bounds = sorted(bounds, key = itemgetter(1)) #if start is shared, sort so shortest isoform is first
            elif gene.strand == '-':
                if bounds[0][1] != bounds[1][1]:
                    bounds = sorted(bounds, key = itemgetter(1), reverse = True)  #sort isoform start/stop by stop coord, highest to lowest
                elif bounds[0][1] == bounds[1][1]:
                    bounds = sorted(bounds, key = itemgetter(0), reverse = True)

        elif event_type == 'ALE':
            if gene.strand == '+':
                if bounds[0][1] != bounds [1][1]:
                    bounds = sorted(bounds, key = itemgetter(1), reverse = True)
                elif bounds[0][1] == bounds [1][1]:
                    bounds = sorted(bounds, key = itemgetter(0), reverse = True)
            elif gene.strand == '-':
                if bounds[0][0] != bounds[1][0]:
                    bounds = sorted(bounds, key = itemgetter(0))
                elif bounds[0][0] == bounds[1][0]:
                    bounds = sorted(bounds, key = itemgetter(1))

    if isoform_to_pick == 'distal':
        bounds_to_pick = bounds[0] #take the first isoform
        start = bounds_to_pick[0]
        stop = bounds_to_pick[1]
        for isoform in db.children(gene, featuretype = 'mRNA'):
            exonlengths = sum(len(i) for i in db.children(isoform, featuretype = 'exon'))
            isoformdata = [isoform.chrom, isoform.start, isoform.stop, isoform.strand, exonlengths]
            if start == isoform.start and stop == isoform.stop and isoformdata not in pickedisoforms:
                pickedisoforms.append(isoformdata)
                print isoform
                for exon in db.children(isoform, featuretype = 'exon'):
                    print exon

    elif isoform_to_pick == 'proximal':
        bounds_to_pick = bounds[1] #take the second isoform
        start = bounds_to_pick[0]
        stop = bounds_to_pick[1]
        for isoform in db.children(gene, featuretype = 'mRNA'):
            exonlengths = sum(len(i) for i in db.children(isoform, featuretype = 'exon'))
            isoformdata = [isoform.chrom, isoform.start, isoform.stop, isoform.strand, exonlengths]
            if start == isoform.start and stop == isoform.stop and isoformdata not in pickedisoforms:
                pickedisoforms.append(isoformdata)
                print isoform
                for exon in db.children(isoform, featuretype = 'exon'):
                    print exon
                



    
os.remove(db_fn)
            

        
