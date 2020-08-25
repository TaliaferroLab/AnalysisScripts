import os
import sys
import argparse
from Bio import SeqIO
import numpy as np
from spats import getkmerpositions
from spats import mappingfilter


def getreactivities(kmerpos, filteredRNAs, reactivities, k):
    #Get reactivities at kmers of interest in RNAs that passed coverage filter
    reactivitydict = {} #{RNAid : {position : reactivity}}
    treatedmodsdict = {} #{RNAid : {position : treated_mods}}
    reactivityfh = open(reactivities, 'r')
    for line in reactivityfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence' and line[6] != '-': #skip header and position 0 of every sequence
            RNAid = line[0]
            if RNAid not in reactivitydict:
                reactivitydict[RNAid] = {}
            if RNAid not in treatedmodsdict:
                treatedmodsdict[RNAid] = {}
            position = int(line[2])
            treated_mods = int(line[4])
            #theta = float(line[7])
            #reactivitydict[RNAid][position] = theta
            bp_prob = float(line[10])
            reactivitydict[RNAid][position] = bp_prob
            treatedmodsdict[RNAid][position] = treated_mods
    reactivityfh.close()

    kmerreactivitydict = {} #{RNAid : {kmerstart : {pos1 : reactivity, pos2 : reactivity, posn : reactivity}}
    for RNA in reactivitydict:
        if RNA in filteredRNAs and RNA in kmerpos: #if it passes both read count filter and it has at least one kmer instance in it
            kmerreactivitydict[RNA] = {}
            for kmerinstance in kmerpos[RNA]:
                ntreactivities = []
                kmertreatedmods = 0
                kmerstart = kmerinstance
                kmerend = kmerinstance + k
                kmerreactivitydict[RNA][kmerstart] = {}
                posinkmer = 1
                for nt in range(kmerstart, kmerend):
                    treated_mods = treatedmodsdict[RNA][nt]
                    kmertreatedmods += treated_mods
                for nt in range(kmerstart, kmerend):
                    reactivity = reactivitydict[RNA][nt]
                    ntreactivities.append(reactivity)
                    kmerreactivitydict[RNA][kmerstart]['pos{0}'.format(posinkmer)] = reactivity
                    posinkmer +=1

    return kmerreactivitydict



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reactivities', type = str, help = 'SHAPE reactivities file containing RNAstructure-calculated bp probabilities.')
    parser.add_argument('--kmers', type = str, help = 'Comma-separated list of kmers (Us not Ts).  However, for this purpose it probably makes more sense to consider one kmer at a time.')
    parser.add_argument('--RNAs', type = str, help = 'RNA sequences in fasta format.')
    parser.add_argument('--mappingfilter', type = int, help = 'RNA stop count filter. Only consider RNAs with at least this many stops in the + channel.')
    parser.add_argument('--rvalues', type = str, help = 'R values from BNS.')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    kmers = args.kmers.split(',')
    kmerpositions = getkmerpositions(kmers, args.RNAs)
    filteredRNAs = mappingfilter(args.reactivities, args.mappingfilter)
    k = int(len(kmers[0]))
    kmerreactivitydict = getreactivities(kmerpositions, filteredRNAs, args.reactivities, k)

    rfh = open(args.rvalues, 'r')
    rs = {}
    for line in rfh:
        line = line.strip().split('\t')
        if line[0] == 'species':
            continue
        RNA = line[1]
        R = float(line[2])
        if R < 2:
            rs[RNA] = 'unbound'
        elif R >= 2:
            rs[RNA] = 'bound'
    rfh.close()

    positions = []
    for i in range(k):
        positions.append('pos{0}'.format(i+1))
    outfh = open(args.outfile, 'w')
    outfh.write(('\t').join(['RNA', 'kmerstart', 'age', 'location', 'BNS'] + positions) + '\n')
    for RNA in kmerreactivitydict:
        age = RNA.split('|')[0]
        location = RNA.split('|')[-1]
        for kmerpos in sorted(kmerreactivitydict[RNA]):
            posreactivities = []
            for position in positions:
                posreactivities.append(str(kmerreactivitydict[RNA][kmerpos][position]))
            outfh.write(('\t').join([RNA, str(kmerpos), age, location, rs[RNA]] + posreactivities) +'\n')

    outfh.close()
