#Reorders MISO Gffs by position of their start coord (if on + strand) or stop coord (if on - strand).  Reorders so that most gene-distal isoform is listed first, followed by all of its children exons.

#Assumes each event (gene) only has two isoforms (mRNAs)!!!

#Usage: python ReorderGff3v2.0.py <gff_file> <event_type (AFE or ALE)> > output

import sys
import os
import gffutils
from operator import itemgetter

gff_fn = sys.argv[1]
event_type = sys.argv[2]
db_fn = os.path.basename(gff_fn) + '.db'

if os.path.isfile(db_fn) == False : #if database doesn't exist, create it
    gffutils.create_db(gff_fn, db_fn)

db = gffutils.FeatureDB(db_fn)

genes= db.features_of_type('gene')


for gene in genes:
    print gene #print gene line
    bounds = []
    for isoform in db.children(gene, featuretype = 'mRNA'):

        if gene.strand == '+':
            bounds.append([isoform.start, isoform.stop])
            
        elif gene.strand == '-':
            bounds.append([isoform.start, isoform.stop])

    if len(bounds) == 1: #if there is only one isoform in event skip it
        continue

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


    for bound in bounds:
        if gene.strand == '+':
            for isoform in db.children(gene, featuretype = 'mRNA'):
                start = bound[0]
                stop = bound[1]
                if start == isoform.start and stop == isoform.stop:
                    print isoform
                    for exon in db.children(isoform, featuretype = 'exon'):
                        print exon

        elif gene.strand == '-':
            for isoform in db.children(gene, featuretype = 'mRNA'):
                start = bound[0]
                stop = bound[1]
                if start == isoform.start and stop == isoform.stop:
                    print isoform
                    for exon in db.children(isoform, featuretype = 'exon'):
                        print exon
    

os.remove(db_fn)

