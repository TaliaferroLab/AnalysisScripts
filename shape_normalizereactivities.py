#This script normalizes reactivity (theta) values from a reactivities.out file produced by spats v 0.8.0.
#It does so following the method outlined in Lucks et al PNAS (2011).
#The top 2% of thetas are excluded.
#Then the 3-10th percentiles of thetas are averaged and all theta values are then normalized by this value.
#This script will also filter out any RNAs that have less than <filter> number of stops in 
#the + channel.  This is to remove RNAs for which there were very few reads.

#The reactivities.out file from spats supplied to this script must be sorted by RNA name (field 1)
#and nucleotide sequence (field 3).  This is the default output.

#Also, when mode == fastasplit, this will split a fasta file of many sequences in many individual fasta
#files of one sequence each. This is useful in preparation for submitting them to the partition function
#in RNAstructure.  If mode == SHAPEfiles, puts SHAPE reactivity files for use with RNAstructure into a separate
#directory.  These fasta and SHAPEfile directories are then ready for use with shape_rnastructure.py

import argparse
from numpy import mean
import os
from Bio import SeqIO
import sys

def filterRNAbystops(reactivities, minstops):
    print 'Filtering RNAs by the number of stops in the + channel...'
    minstops = int(minstops)
    reactivitiesfh = open(reactivities, 'r')
    stopsdict = {} # {RNA : number of cumulative stops in + channel at all nucleotide positions}
    filteredRNAs = [] #list of RNAs that pass the stops filter
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence': #skip header
            sequence = line[0]
            treatedmods = int(line[4])
            if sequence in stopsdict:
                stopsdict[sequence] += treatedmods
            elif sequence not in stopsdict:
                stopsdict[sequence] = treatedmods

    for RNA in stopsdict:
        if stopsdict[RNA] >= minstops:
            filteredRNAs.append(RNA)

    print '{0} of {1} RNAs have at least {2} stops in the treated channel.'.format(len(filteredRNAs), len(stopsdict), minstops)
    reactivitiesfh.close()

    return filteredRNAs

def getThetanormfactor(reactivities):
    print 'Normalizing theta values according to the 2%/8% rule...'
    reactivitiesfh = open(reactivities, 'r')
    reactivities = [] #list of all reactivity values, unnormalized
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        #skip header and the 5' most nucleotide of all RNAs, which has theta = '-' and any nucleotide after 108, which never has any stops in either channel
        if line[0] != 'sequence' and int(line[2]) > 0 and int(line[2]) < 109:
            theta = float(line[7])
            reactivities.append(theta)

    reactivities = sorted(reactivities, reverse=True)
    top2percentile = int(len(reactivities) * 0.02)
    top10percentile = int(len(reactivities) * 0.1)
    print 'There were {0} reactivities. The top 2% are therefore the first {1} and the top 10% are the first {2}.'.format(len(reactivities), top2percentile, top10percentile)
    #Calculate average of 3-10 percentile
    normfactor = float(mean(reactivities[top2percentile + 1 : top10percentile]))
    print 'After throwing out the top 2%, the average reactivity of the next 8% is {0}.'.format(normfactor)

    reactivitiesfh.close()

    return normfactor

def writeNormalizedreactivities(reactivities, normfactor, outfile):
    reactivitiesfh = open(reactivities, 'r')
    outfh = open(outfile, 'w')
    outfh.write(('\t').join(['sequence','rt_start','five_prime_offset','nucleotide','treated_mods','untreated_mods','beta','theta','normalized_theta','c']) + '\n')
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] == 'sequence':
            continue
        if line[7] == '-': #if this is the 5' most nucleotide and therefore theta = '-'
            outfh.write(('\t').join([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], '-', line[8]]) + '\n')
            continue
        outfh.write(('\t').join([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], str(float(line[7]) / normfactor), line[8]]) + '\n')

    outfh.close()


def writeSHAPEfiles(reactivities, normfactor, filteredRNAs, outdir):
    reactivitiesfh = open(reactivities, 'r')
    #See if outdir exists.  If not, make it.
    if os.path.isdir(os.path.abspath(outdir)) == False:
        os.mkdir(os.path.abspath(outdir))
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence':
            RNAname = line[0]
            nucleotide = int(line[2])
            if nucleotide == 0: #if this is the first time we see this RNA
                outfh = open(os.path.join(os.path.abspath(outdir), RNAname + '.SHAPE'), 'w')
            elif nucleotide > 0 and nucleotide < 133: #any nucleotide after 133 has a reactivity of 0 due to RT primer binding
                theta = float(line[7])
                normalizedtheta = theta/normfactor
                outfh.write(str(nucleotide) + '\t' + str(normalizedtheta) + '\n')
            elif nucleotide == 133:
                outfh.close()
            elif nucleotide > 133:
                continue

def splitfasta(fasta, outdir):
    #Split a fasta file of many sequences into many small files of one sequence each.
    #See if outdir exists. If not, make it.
    if os.path.isdir(os.path.abspath(outdir)) == False:
        os.mkdir(os.path.abspath(outdir))

    for record in SeqIO.parse(fasta, 'fasta'):
        outfh = open(os.path.join(os.path.abspath(outdir), record.id + '.fasta'), 'w')
        outfh.write('>' + str(record.id) + '\n')
        outfh.write(str(record.seq))
        outfh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reactivities', type = str, help = 'reactivities.out from spats output.')
    parser.add_argument('--outfile', type = str, help = 'Required if writing normalized reactivities.out.')
    parser.add_argument('--minstops', type = str, help = 'Required if writing SHAPE files. Minimum number of stops required in the + channel for an RNA to be considered.')
    parser.add_argument('--fasta', type = str, help = 'Required if using fastasplit mode. A fasta file of many sequences to be split to many small files of one sequence each.')
    parser.add_argument('--outdir', type = str, help = 'Required if writing SHAPE reactivity files for RNAstructure OR splitting fasta file into many files with a single sequence.')
    parser.add_argument('--mode', type = str, required = True, choices = ['normalizetheta', 'SHAPEfiles', 'splitfasta'], 
                        help = 'Choose normalizetheta to output a reactivities.out file with normalized theta values. Choose SHAPEfiles to output many SHAPE reactivity files for use with RNA structure. Choose fastasplit to split a fasta file of many sequences into many files of one sequence each.')
    args = parser.parse_args()

    if args.mode == 'normalizetheta':
        if args.reactivities == None:
            print 'Error: must provide reactivities.out from spats output.'
            sys.exit()
        normfactor = getThetanormfactor(args.reactivities)
        writeNormalizedreactivities(args.reactivities, normfactor, args.outfile)

    elif args.mode == 'SHAPEfiles':
        if args.reactivities == None:
            print 'Error: must provide reactivities.out from spats output.'
            sys.exit()
        filteredRNAs = filterRNAbystops(args.reactivities, args.minstops)
        normfactor = getThetanormfactor(args.reactivities)
        writeSHAPEfiles(args.reactivities, normfactor, filteredRNAs, args.outdir)

    elif args.mode == 'splitfasta':
        if args.fasta == None or args.outdir == None:
            print 'Error: must provide a fasta to split and a directory to put split files in.'
            sys.exit()
        splitfasta(args.fasta, args.outdir)
    

