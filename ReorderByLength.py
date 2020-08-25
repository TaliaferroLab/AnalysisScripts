#reorders MISO gff files based on mRNA lengths.  mRNA lengths are the combined lengths of their children exons.  Exons themselves are not sorted

#Usage: python ReorderByLength.py [infile.gff3] > outfile.gff3

import sys
import os
import gffutils
import operator

gff_fn = sys.argv[1]
db_fn = os.path.basename(gff_fn) + '.db'

if os.path.isfile(db_fn) == False : #if database doesn't exist, create it
    gffutils.create_db(gff_fn, db_fn)

db = gffutils.FeatureDB(db_fn)

genes= db.features_of_type('gene')

for gene in genes:
    print gene
    lengths = {}
    for isoform in db.children(gene, featuretype='mRNA'):
        lengths[isoform.id] = sum(len(i) for i in db.children(isoform, featuretype='exon')) # get isoform lengths
        sortedlengths = sorted(lengths.iteritems(), key= operator.itemgetter(1), reverse=True) #sorted lengths longest to shortest...sorted lengths is a list of tuples of isoform IDs and their lengths
        
    for sortedlength in sortedlengths:
        isoformid = sortedlength[0]
        for isoform in db.children(gene, featuretype='mRNA'):
            if isoform.id == isoformid: #if the first item in the tuple matches any isoform.id
                print isoform

    for exon in db.children(gene, featuretype='exon'):
        print exon #exons are not sorted
        

