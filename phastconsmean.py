#Take a gff (of, say, 3' UTRs) and determine the mean conservation of each base in that 3' UTR.  Each UTR will thus get one "score".
#Gff's are 0-based and phastcons beds are also 0-based

#Usage: python phastconsmean.py <gfffileof3'UTRs> <phastconsbedfile> <UTRtype> <outfile>

import sys
import os
import gffutils
import tabix
import pysam

def phastconsmean(gff, phastconsbed):
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)

    phastconsmeans = {} # {UTR_ID : mean phastcons score}
    interrogated_3UTRs = 0
    interrogated_bps = 0
    UTRs_with_scores = 0
    bps_with_scores = 0
    phastconstabix = pysam.Tabixfile(phastconsbed)
    three_prime_UTRs = db.features_of_type('UTR3')

    for UTR in three_prime_UTRs:
        phastconsscores = []
        interrogated_3UTRs +=1
        if interrogated_3UTRs % 100 == 0:
            sys.stderr.write('Interrogating UTR {0}\n'.format(interrogated_3UTRs))
        for exon in db.children(UTR, featuretype = 'exon', level = 1):
            #for bed in phastconstabix.fetch(str(UTR.chrom), exon.start, exon.end, parser = pysam.asBed()):
            for bed in phastconstabix.fetch(str(UTR.chrom), exon.start, exon.end, parser = pysam.asTuple()):
                score = bed[4]
                #phastconsscores.append(float(bed.name))
                phastconsscores.append(float(score))

        if len(phastconsscores) > 0:
            UTRs_with_scores +=1
            phastconsmeans[UTR.id] = (sum(phastconsscores) / float(len(phastconsscores)))
            
        UTRlength = int(UTR.stop - UTR.start)
        scored_bps = len(phastconsscores)
        interrogated_bps += UTRlength
        bps_with_scores += scored_bps

    percent_bp_with_scores = (bps_with_scores / float(interrogated_bps))*100

    sys.stderr.write('Searched {0} UTRs and got PhastCons scores for {1} of them. Overall, scores were obtained for {2}% of all interrogated bps.\n'.format(interrogated_3UTRs, UTRs_with_scores, percent_bp_with_scores))

    os.remove(db_fn)

    return phastconsmeans

if __name__ == '__main__':
    outfh = open(sys.argv[4], 'w')
    outfh.write('UTR' + '\t' + 'PhastCons_mean' + '\t' + 'Class' + '\n')
    phastconsmeans = phastconsmean(sys.argv[1], sys.argv[2])
    for UTR in phastconsmeans:
        outfh.write(str(UTR) + '\t' + str(phastconsmeans[UTR]) + '\t' + str(sys.argv[3]) + '\n')

    outfh.close()

            
