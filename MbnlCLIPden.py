#Takes a gff of UTRs and a bed file of CLIP clusters and intersects them.  Returns density of CLIP clusters in each UTR.

#Usage: python MbnlCLIPden.py <gff_file> <clusters_bed_file> <outfile.txt>

import sys
import os
import tabix
import pysam
import gffutils
import argparse

def getFPKMs(FPKMfile):
    FPKMdict = {} #{event : FPKM of gene}
    infh = open(FPKMfile, 'r')
    for line in infh:
        line = line.strip().split('\t')
        event = line[0]
        FPKM = float(line[2])
        FPKMdict[event] = FPKM
    infh.close()
    return FPKMdict

def getDeltaDeltaPsis(DeltaDeltaPsis):
    DeltaDeltaPsidict = {} # {event : deltadeltapsi}
    infh = open(DeltaDeltaPsis, 'r')
    for line in infh:
        line = line.strip().split('\t')
        event = line[0]
        deltadeltapsi = str(line[8])
        DeltaDeltaPsidict[event] = deltadeltapsi
    infh.close()
    return DeltaDeltaPsidict

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
            print 'Warning: GFF region has start {0} > end {1}.  Skipping...'.format(UTR.start, UTR.stop)
            continue
        try:
            cluster_density = float((UTR_intersecting_clusters) / (UTRlength))
        except ZeroDivisionError:
            print 'WARNING: UTR {0} has length 0.'.format(UTR.id)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type = str, help = 'Output file.', required = True)
    parser.add_argument('--gff', type = str, help = 'Gff of regions to interrogate.', required = True)
    parser.add_argument('--clusters', type = str, help = 'Tabix indexed file of clip clusters.', required = True)
    parser.add_argument('--FPKMs', type = str, help = 'File of FPKMs for genes containing events.  Tab separated. Event, Gene, FPKM.', required = False)
    parser.add_argument('--deltadeltapsis', type = str, help = 'Table of delta delta psis.', required = False)
    args = parser.parse_args()
    
    if args.FPKMs == None and args.deltadeltapsis == None:
        outfh = open(args.output, 'w')
        outfh.write('UTR' + '\t' + 'clusters' + '\t' + 'UTRLength' + '\t' + 'density' + '\n')
        cluster_densities = MbnlCLIPden(args.gff, args.clusters)
        for UTR in cluster_densities:
            clusters = str(cluster_densities[UTR][0])
            UTRLength = str(cluster_densities[UTR][1])
            density = str(cluster_densities[UTR][2])
            ###
            #outfh.write(UTR + '\t' + clusters + '\t' + UTRLength + '\t' + density + '\n')
            #THIS HAS BEEN CHANGED TO GIVE THE NAME OF THE EVENT, NOT THE ISOFORM
            ###
            outfh.write(UTR[:-2] + '\t' + clusters + '\t' + UTRLength + '\t' + density + '\n')
        outfh.close()

    elif args.FPKMs and args.deltadeltapsis:
        outfh = open(args.output, 'w')
        outfh.write('UTR' + '\t' + 'clusters' + '\t' + 'UTRLength' + '\t' + 'density' + '\t' + 'FPKM' + '\t' + 'FPKMnormalizeddensity' + '\t' + 'deltadeltapsi' + '\t' + 'densitybin' + '\n')
        cluster_densities = MbnlCLIPden(args.gff, args.clusters)
        FPKMs = getFPKMs(args.FPKMs)
        deltadeltapsis = getDeltaDeltaPsis(args.deltadeltapsis)
        for UTR in cluster_densities:
            if UTR[:-2] in FPKMs:
                densitybin = ''
                clusters = str(cluster_densities[UTR][0])
                UTRLength = str(cluster_densities[UTR][1])
                density = str(cluster_densities[UTR][2])
                FPKM = float(FPKMs[UTR[:-2]])
                FPKMnormalizeddensity = float(density) / FPKM
                deltadeltapsi = deltadeltapsis[UTR[:-2]]
                if float(density) == 0:
                    densitybin = 'none'
                elif float(density) >0 and float(density) < 0.002:
                    densitybin = 'low'
                elif float(density) >= 0.002 and float(density) < 0.006:
                    densitybin = 'medium'
                elif float(density) >= 0.006:
                    densitybin = 'high'
                outfh.write(UTR + '\t' + clusters + '\t' + UTRLength + '\t' + density + '\t' + str(FPKM) + '\t' + str(FPKMnormalizeddensity) + '\t' + deltadeltapsi + '\t' + densitybin + '\n')

        outfh.close()

    elif args.FPKMs == None and args.deltadeltapsis:
        outfh = open(args.output, 'w')
        outfh.write('UTR' + '\t' + 'clusters' + '\t' + 'UTRLength' + '\t' + 'density' + '\t' + 'deltadeltapsi' + '\n')
        cluster_densities = MbnlCLIPden(args.gff, args.clusters)
        deltadeltapsis = getDeltaDeltaPsis(args.deltadeltapsis)
        print cluster_densities
        for UTR in cluster_densities:
            clusters = str(cluster_densities[UTR][0])
            UTRLength = str(cluster_densities[UTR][1])
            density = str(cluster_densities[UTR][2])
            if UTR[:-2] in deltadeltapsis:
                deltadeltapsi = deltadeltapsis[UTR[:-2]]
                outfh.write(UTR + '\t' + clusters + '\t' + UTRLength + '\t' + density + '\t' + deltadeltapsi + '\n')

        outfh.close()
