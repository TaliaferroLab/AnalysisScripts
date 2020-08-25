#Takes a directory of PSSMs (usually ~/Documents/MIT/Localization/mousekmertoprotein/pwms_all_motifs), a file correlating PSSMids with protein ids
#(usually ~/Documents/MIT/Localization/mousekmertoprotein/RBP_Information.txt),
#and k.  Slides kmer 1 nt at a time along the PSSM and takes the highest scoring window as the score for that 
#kmer and PSSM. Makes all possible kmers of length k and then get scores for every kmer against every protein.  To be considered
#a match for that protein, the kmer must be in the top X percentage of matches for all kmers (set by "threshhold") and it must meet
#that filter in both biological replicates for the protein ('setA' and 'setB').

import os
import sys
import argparse
import operator
import itertools

def motiftoprotein(masterfile):
    #Convert between motif names (e.g. RNCMPT00106) and protein names (e.g. SRSF1)
    motiftoprotein = {} #{motifid : [proteinid]}

    masterfh = open(masterfile, 'r')
    for line in masterfh:
        line = line.strip().split('\t')
        motifid = line[3]
        proteinid = line[6]
        species = line[7]
        motiftoprotein[motifid] = [proteinid, species]
    masterfh.close()

    return motiftoprotein

def getPSSMs(directory):
    PSSMs = {} # {protein : [{A : 0.1, C : 0.1, G : 0.8, U : 0.8}, {A : 0.2, C : 0.3, G : 0.5, U : 0.5}], protein2 : ...}
    for file in os.listdir(directory):
        if os.path.basename(file).endswith('.txt'): #PSSMs must be the only thing in the directory that end with .txt!!!
            filename = os.path.basename(file).split('.txt')[0]
            PSSMs[filename] = []
            
            fh = open(os.path.join(directory, file), 'r')
            for line in fh:
                line = line.strip().split('\t')
                if 'Pos' not in line: #skip header
                    PSSMs[filename].append({'A' : float(line[1]), 'C' : float(line[2]), 'G' : float(line[3]), 'U' : float(line[4])})

            fh.close()
    
    return PSSMs

def makeallkmers(k):
    k = int(k)
    bases = ['A','U','G','C']
    print 'Making all possible {0}mers.'.format(k)
    allkmers = [''.join(x) for x in itertools.product(bases, repeat=k)]
    print 'Done!'

    return allkmers

def getkmerscore(kmer, PSSM):
    kmerlen = len(kmer)
    startpos = 0
    kmerscores = []
    PSSMlength = len(PSSM.values()[0]) #number of positions in PSSM
    
    #Some PSSM files are empty.  Dunno why.  Also, some PSSMs are shorter than the length of the kmer
    if PSSMlength == 0 or PSSMlength < kmerlen:
        return None
    while startpos <= (PSSMlength - kmerlen):
        ntscores = [] #reset list
        for freqs in PSSM.values(): #slide kmer over 1 by 1, fitting the kmer in the PSSM as many times as you can
            for i in range(0, kmerlen):
                ntscores.append(freqs[i + startpos][kmer[i]])
            kmerscore = reduce(operator.mul, ntscores, 1) #kmerscore is all ntscores multiplied together
            kmerscores.append(kmerscore)
            
            startpos +=1
    
    return max(kmerscores) #return the highest score

def getpassingkmers(PSSMs, allkmers, threshhold, motiftoprotein):
    #threshhold is fraction between 0 and 1. Lower threshhold is more stringent.
    threshhold = float(threshhold)
    passingkmers = {} #{motif : [passing kmers]}
    pwm = 0
    for motif in PSSMs:
        #Skip any motif that doesn't have a protein match
        if motif not in motiftoprotein:
            continue
        proteinname = motiftoprotein[motif][0]
        species = motiftoprotein[motif][1]
        passingkmers[proteinname] = []
        PSSM = {motif : PSSMs[motif]}
        scores = {} # {kmer : score}
        for kmer in allkmers:
            score = getkmerscore(kmer, PSSM)
            if score: #Ignore anything that returned None.
                scores[kmer] = score

        no_of_passes = len(allkmers) * threshhold #this is likely a float
        sortedscores = sorted(scores.items(), key = operator.itemgetter(1), reverse = True)
        #Scores are now a list of tuples (kmer, score), sorted in descending order
        for ind, tup in enumerate(sortedscores):
            #Because index is 0 based, gotta add one
            if ind + 1 <= no_of_passes and (species == 'Mus_musculus'):
                passingkmers[proteinname].append(tup[0])

    return passingkmers

def getpassingproteins(passingkmers, allkmers):
    #passingkmers is protein centric. We need something that is kmer centric instead. Just reorganize the dictionary.
    passingproteins = {} #{kmer : [list of passingproteins]}
    for protein in passingkmers:
        kmers = passingkmers[protein]
        for kmer in kmers:
            if passingproteins.has_key(kmer) == True:
                passingproteins[kmer].append(protein)
            elif passingproteins.has_key(kmer) == False:
                passingproteins[kmer] = [protein]

    #Any kmer that doesn't match to any protein still needs to be in this dictionary
    for kmer in allkmers:
        if kmer not in passingproteins:
            passingproteins[kmer] = ['None']

    return passingproteins
    


if __name__ == '__main__':
    motiftoprotein = motiftoprotein(sys.argv[1])
    PSSMs = getPSSMs(sys.argv[2])
    allkmers = makeallkmers(sys.argv[3])
    passingkmers = getpassingkmers(PSSMs, allkmers, sys.argv[4], motiftoprotein)
    getpassingproteins(passingkmers, allkmers)
