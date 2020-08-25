import sys
import os
import tabix
import pysam
import gffutils
import argparse
import random
import numpy
from scipy import stats


def MbnlCLIPden(gff, clusters):
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, force=True)

    db = gffutils.FeatureDB(db_fn)

    cluster_densities = {} # {UTR_ID : [number of intersecting clusters, UTR length, cluster density in that UTR]}
    UTRlengths = []
    number_of_UTRs = len(list(db.features_of_type('3\'UTR')))
    print number_of_UTRs
    three_prime_UTRs = db.features_of_type('3\'UTR')
    clustertabix = pysam.Tabixfile(clusters)
    total_intersecting_clusters = 0

    for UTR in three_prime_UTRs:
        UTR_intersecting_clusters = 0
        UTRlength = float(UTR.stop - UTR.start)
        UTRlengths.append(UTRlength)
        try:
            for bed in clustertabix.fetch(str(UTR.chrom), UTR.start, UTR.stop, parser = pysam.asBed()):
                if bed.strand == UTR.strand:
                    '''
                    print UTRlength, UTR.id, bed.start, bed.end
                    '''
                    UTR_intersecting_clusters +=1
                    total_intersecting_clusters +=1
        except ValueError:
            #print 'Warning: GFF region has start {0} > end {1}.  Skipping...'.format(UTR.start, UTR.stop)
            continue
        try:
            cluster_density = float((UTR_intersecting_clusters) / (UTRlength))
        except ZeroDivisionError:
            #print 'WARNING: UTR {0} has length 0.'.format(UTR.id)
            continue
        cluster_densities[str(UTR.id)] = [UTR_intersecting_clusters, UTRlength, cluster_density]

    total_sequence_length = sum(UTRlengths)

    sys.stderr.write('Searched {0} nt from {1} UTRs and found {2} clusters within them. The overall cluster density is {3} clusters per nt \n'.format(int(total_sequence_length), number_of_UTRs, total_intersecting_clusters, float(float(total_intersecting_clusters) / float(total_sequence_length))))

    os.remove(db_fn)
    if os.path.isfile(db_fn + '-shm'):
        os.remove(db_fn + '-shm')
    if os.path.isfile(db_fn + '-wal'):
        os.remove(db_fn + '-wal')

    return cluster_densities


def samplevalues(distaldensities, proximaldensities, number_to_sample, times_to_sample, threshold):
    #Given a set of clip cluster densities (distal and proximal) defined by the above function (in this case likely from
    #a large set of "unaffected" UTRs, sample it some number of times, recording the mean clip density
    #every time. Count what fraction of time this mean exceeds a "threshold", which can be the mean
    #clip cluster density for "affected" UTRs.

    #densitydict format = # {UTR_ID : [number of intersecting clusters, UTR length, cluster density in that UTR]}


    ALEids = []
    commonids = []
    relativedensities = []
    greater_than_threshold = 0

    #commonids have to be present in distaldensites and proximaldensities
    for UTR_ID in distaldensities:
        #UTR_IDs are XXXX.A or XXXX.B, depending on which isoform they came from. The name of the parent ALE is therefore
        #XXXX. Get this name.
        ALEids.append(('.').join(UTR_ID.split('.')[:-1]))
    for UTR_ID in proximaldensities:
        if ('.').join(UTR_ID.split('.')[:-1]) in ALEids:
            commonids.append(('.').join(UTR_ID.split('.')[:-1]))

    print 'From {0} distal UTRs and {1} proximal UTRs, found {2} in common that make up ALEs.'.format(len(distaldensities), len(proximaldensities), len(commonids))

    for i in range(int(times_to_sample)):
        randomids = []
        randomids = random.sample(commonids, int(number_to_sample))
        distalclusters = 0
        distallength = 0
        proximalclusters = 0
        proximallength = 0
        for ALEID in randomids:
            for UTRID in distaldensities:
                if ALEID in UTRID:
                    distalclusters += distaldensities[UTRID][0]
                    distallength += distaldensities[UTRID][1]

            for UTRID in proximaldensities:
                if ALEID in UTRID:
                    proximalclusters += proximaldensities[UTRID][0]
                    proximallength += proximaldensities[UTRID][1]

        distaldensity = distalclusters / float(distallength)
        proximaldensity = proximalclusters / float(proximallength)
        if proximaldensity == 0:
            continue
        relativedensity = distaldensity / proximaldensity
        relativedensities.append(relativedensity)
        if i % 100 == 0:
            print 'Random sample {0}...'.format(i)

    for density in relativedensities:
        if density > threshold:
            greater_than_threshold +=1

    fractiongreater = greater_than_threshold / float(len(relativedensities))
                    

    print 'Of {0} random samplings of {1} UTRs each, {2} ({3}) had mean relative densities greater than the threshold.'.format(times_to_sample, number_to_sample, 
                                                                                                                               greater_than_threshold, fractiongreater)
    print 'The mean relative density was {0} with a SD of {1} and a SE of {2}.'.format(numpy.mean(relativedensities), numpy.std(relativedensities), stats.sem(relativedensities))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--distalgff', type = str, help = 'Gff of distal regions to interrogate.', required = True)
    parser.add_argument('--proximalgff', type = str, help = 'Gff of proximal regions to interrogate.', required = True)
    parser.add_argument('--clusters', type = str, help = 'Tabix indexed file of clip clusters.', required = True)
    parser.add_argument('--numbertosample', type = int, help = 'Number of relative densities to sample in each sampling.', required = True)
    parser.add_argument('--timestosample', type = int, help = 'Number of times to sample relative densities.', required = True)
    parser.add_argument('--threshold', type = float, help = 'Threshold relative density to compare against.', required = True)
    args = parser.parse_args()

    distaldensities = MbnlCLIPden(args.distalgff, args.clusters)
    proximaldensities = MbnlCLIPden(args.proximalgff, args.clusters)

    samplevalues(distaldensities, proximaldensities, args.numbertosample, args.timestosample, args.threshold)
