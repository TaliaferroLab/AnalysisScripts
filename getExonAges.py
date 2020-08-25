#Takes a gff of exons you are interested in and a table of exon ages.  This table was given to me by Jason
#and has the following format:
#
#10:100057171:100057290:+	300	110 
#where the 2nd field is the genomic age and the 3rd field is the splicing age.
#Returns a file of exons in your file that have ANY OVERLAP with an exon in the table of exon ages.
#This strategy was used because (per Jason) the boundaries of the exons in the ages table
#are not well defined.

#Usage python getExonAges.py --help

import sys
import argparse
import gffutils
import os

def getExonAges(exons, exonagestable):
    allexonages = {} #{(chrm, start, stop, strand) : [genomicage, splicingage]}
    selectedexonages = {} #{ID : [genomicage, splicingage]}
    allexonagesfh = open(exonagestable, 'r')
    queriedexons = 0
    for line in allexonagesfh:
        line = line.strip().split('\t')
        if 'GenomicAge' not in line: #skip header
            chrm = 'chr' + line[0].split(':')[0]
            start = int(line[0].split(':')[1])
            stop = int(line[0].split(':')[2])
            strand = line[0].split(':')[3]
            genomicage = line[1]
            splicingage = line[2]
            allexonages[(chrm, start, stop, strand)] = [genomicage, splicingage]

    allexonagesfh.close()

    gff_fn = exons
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn): #if there's already a database there, remove it
        os.remove(db_fn)

    gffutils.create_db(gff_fn, db_fn)
    db = gffutils.FeatureDB(db_fn)

    features = db.features_of_type('mRNA') #CHANGE THIS TO MATCH 3RD FIELD OF GFF
    exonstoquery =  sum(1 for _ in features)
    features = db.features_of_type('mRNA') #remake generator

    for feature in features:
        queriedexons +=1
        if queriedexons % 10 == 0:
            print 'Querying exon {0} of {1}...\n'.format(queriedexons, exonstoquery)
        featurerange = range(feature.start, feature.stop)
        for exon in allexonages:
            genomicage = allexonages[exon][0]
            splicingage = allexonages[exon][1]
            exonrange = range(exon[1], exon[2])
            overlap = list(set(featurerange) & set(exonrange))
            if feature.chrom == exon[0] and feature.strand == exon[3] and len(overlap) > 0:
                selectedexonages[feature.id] = [genomicage, splicingage]

    os.remove(db_fn)

    return selectedexonages

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exons', type = str, help = 'Gff of exons to consider. Usually miso annotations.')
    parser.add_argument('--exonages', type = str, help = 'Table of exon ages. Format: chrm:start:stop:strand \t genomicage \t splicingage')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    selectedexonages = getExonAges(args.exons, args.exonages)
    outfh = open(args.outfile, 'w')
    outfh.write('Exon' + '\t' + 'GenomicAge' + '\t' + 'SplicingAge' + '\n')
    for exon in selectedexonages:
        genomicage = selectedexonages[exon][0]
        splicingage = selectedexonages[exon][1]
        outfh.write(exon + '\t' + genomicage + '\t' + splicingage + '\n')
    outfh.close()
