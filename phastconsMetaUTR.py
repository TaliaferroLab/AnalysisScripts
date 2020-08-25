#Make metaUTR plots of phastcons values.  Takes as input a gff of UTRs you are interested in and 
#a .bed file of phastcons values (sorted.phastcons.mm9.bed.gz)


import sys
import os
import argparse
import gffutils
import tabix
import pysam
from numpy import mean

def phastconsMetaUTR(gff, phastconsbed, numberofbins):
    binfactor = 100 / float(numberofbins)
    phastconsdict = {} #{binnumber: [phastconsvalue1, phastconsvalue2, ...]}
    averagedict = {} #{binnumber : meanphastconsscore}

    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    phastconstabix = pysam.Tabixfile(phastconsbed)
    three_prime_UTRs = db.features_of_type('3\'UTR')
    interrogated_3UTRs = 0

    for UTR in three_prime_UTRs:
        interrogated_3UTRs +=1
        if interrogated_3UTRs % 100 == 0:
            sys.stderr.write('Interrogating UTR {0}\n'.format(interrogated_3UTRs))

        UTRlength = float(UTR.end - UTR.start)
        phastconsscores = {} #{base : score}
        #Get every bed entry that overlaps with UTR coords
        try: #Occasionally startcoord > endcoord or you get a UTR on a chrm that is not in the phastcons bed file (like chrY)
            for bed in phastconstabix.fetch(str(UTR.chrom), UTR.start, UTR.stop, parser = pysam.asBed()):
                phastconsscores[float(bed.start)] = float(bed.name)
        except ValueError:
            print 'WARNING: something is wrong with UTR {0}. Either start > end or it\'s on a chromosome that is not in the phastcons bed file.'.format(UTR.id)
            continue
                
        if len(phastconsscores) > 0: #if there were any bases in the UTR that had phastcons scores
            for base in phastconsscores:
                position = base - UTR.start #nth base of UTR
                binnumber = round((float(position / float(UTRlength))/binfactor),2) * binfactor
                if phastconsdict.has_key(binnumber) == False and binnumber >= 0 and binnumber <= 1:
                    phastconsdict[binnumber] = [phastconsscores[base]]
                elif phastconsdict.has_key(binnumber) == True and binnumber >= 0 and binnumber <= 1:
                    phastconsdict[binnumber].append(phastconsscores[base])

    #Calculate mean scores at every bin
    for binnumber in phastconsdict:
        averagedict[binnumber] = mean(phastconsdict[binnumber])

    os.remove(db_fn)
    '''
    os.remove((os.path.basename(db_fn) + '-shm'))
    os.remove((os.path.basename(db_fn) + '-wal'))
    '''

    return averagedict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', type = str, help = 'Gff of UTRs to analyze.')
    parser.add_argument('--phastconsbed', type = str, help = 'Phastcons scores in bed format and gzipped. Must also have tabix index in same directory.')
    parser.add_argument('--numberofbins', type = int, help = 'Number of bins in metagene.')
    parser.add_argument('--output', type = str, help = 'Output file.')
    args = parser.parse_args()

    averagedict = phastconsMetaUTR(args.gff, args.phastconsbed, args.numberofbins)

    outfh = open(args.output, 'w')
    outfh.write('bin' + '\t' + 'mean_phastcons_score' + '\n')
    for binnumber in sorted(averagedict):
        outfh.write(str(binnumber) + '\t' + str(averagedict[binnumber]) + '\n')
    outfh.close()
