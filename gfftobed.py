#Converts gff to bed. 
#If this isn't what you want, make small modifications.

#Usage: python gfftobed.py --help

import gffutils
import os
import argparse

def gff2bed(gff, featuretype, outfile):
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'
    outfh = open(outfile, 'w')

    if os.path.isfile(db_fn) == False: #if database doesn't exist, create it
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)

    features = db.features_of_type(featuretype)

    for feature in features:
        #THIS COULD REMOVES STOP CODONS AND LAST 50 NT OF UTR!!!!
        
        if feature.strand == '+':
            outfh.write(feature.chrom + '\t' + str(feature.start-3) + '\t' + str(feature.stop) + '\t' + 'ID=' + feature.id + '\t' + '1000' + '\t' + feature.strand + '\n')
        elif feature.strand == '-':
            outfh.write(feature.chrom + '\t' + str(feature.start) + '\t' + str(feature.stop+3) + '\t' + 'ID=' + feature.id + '\t' + '1000' + '\t' + feature.strand + '\n')

    outfh.close()
    os.remove(db_fn)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', type = str, help = 'gff file to convert')
    parser.add_argument('--featuretype', type = str, help = 'Feature type in gff to convert. Third field of gff.')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    gff2bed(args.gff, args.featuretype, args.outfile)
