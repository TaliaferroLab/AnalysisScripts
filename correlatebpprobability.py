import os
import argparse
from Bio import SeqIO
import numpy as np

def mappingfilter(mappingstats, filt):
    #Only consider RNAs with <filt> stops in the + channel
    print 'Filtering RNAs by stop counts...'
    filt = int(filt)
    RNAcount = 0
    filteredRNAcount = 0
    filteredRNAs = []
    statsfh = open(mappingstats, 'r')
    for line in statsfh:
        line = line.strip().split('\t')
        if line[0] != 'target': #skip header
            RNAcount +=1
            RNAid = line[0]
            treated_stops = int(line[3])
            if treated_stops >= filt:
                filteredRNAs.append(RNAid)
                filteredRNAcount +=1

    statsfh.close()

    print '{0} of {1} RNAs pass stop count filter.'.format(filteredRNAcount, RNAcount)

    return filteredRNAs

def parsebpprobs(bpprobs): #parse bp probability file from nicole (name and probabilities, tab separated)
    bpprobdict = {} #{RNAname : [list_of_bp_probabilities]
    bpfh = open(bpprobs, 'r')
    for line in bpfh:
        probs = []
        line = line.strip().split('\t')
        for idx, prob in enumerate(line):
            if idx == 0:
                RNAname = prob
            elif idx > 0:
                probs.append(float(prob))
        bpprobdict[RNAname] = probs

    bpfh.close()
    return bpprobdict

def getreactivities(reactivities):
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
            position = int(line[1])
            treated_mods = int(line[3])
            theta = float(line[6])
            reactivitydict[RNAid][position] = theta
            treatedmodsdict[RNAid][position] = treated_mods
    reactivityfh.close()

    return reactivitydict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mappingstats', type = str, help = 'mapping_stats.txt from spats output.', required = True)
    parser.add_argument('--mappingfilter', type = int, help = 'RNA stop count filter. Only consider RNAs with at least this many stops in the + channel.', required = True)
    parser.add_argument('--reactivities', type = str, help = 'reactivities.out from spats run.', required = True)
    parser.add_argument('--output', type = str, help = 'Output file.', required = True)
    parser.add_argument('--bpprobs', type = str, help = 'File of in silico caluclated bp probabilities. Name followed by probabilities.  Values are tab separated.')

    args = parser.parse_args()

    filteredRNAs = mappingfilter(args.mappingstats, args.mappingfilter)
    bpprobdict = parsebpprobs(args.bpprobs)
    reactivitydict = getreactivities(args.reactivities)
    outdict = {} # {RNA : [[pos1, bpprob, reactivity], [pos2, bpprob, reactivity]...]}
    for RNA in filteredRNAs:
        probsandreacts = []
        for i in range(1,129): #Don't consider places where the RT primer was binding
            bpprob = bpprobdict[RNA][i-1]
            reactivity = reactivitydict[RNA][i]
            probsandreacts.append([i, bpprob, reactivity])
        outdict[RNA] = probsandreacts

    outfh = open(args.output, 'w')
    outfh.write('RNA' + '\t' + 'position' + '\t' + 'bpprob' + '\t' + 'reactivity' + '\n')
    for RNA in outdict:
        for pos in outdict[RNA]:
            outfh.write(RNA + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\t' + str(pos[2]) + '\n')

    outfh.close()
    
