#Takes a directory of PSSMs (usually ~/Documents/MIT/Localization/kmertoprotein/), a file correlating PSSMids with protein ids,
#and a kmer of interest.  Slides kmer 1 nt at a time along the PSSM and takes the highest scoring window as the score for that 
#kmer and PSSM.  Returns the top 2% scoring protein matches for the kmer.

#python kmertoprotein.py --help

import os
import sys
import argparse
import operator

def getPSSMs(directory):
    PSSMs = {} # {protein : [{A : 0.1, C : 0.1, G : 0.8, U : 0.8}, {A : 0.2, C : 0.3, G : 0.5, U : 0.5}], protein2 : ...}
    for file in os.listdir(directory):
        if 'top10align' in os.path.basename(file):
            filename = os.path.basename(file).split('_')[0] + '_' + os.path.basename(file).split('_')[2]
            PSSMs[filename] = []
            
            fh = open(os.path.join(directory, file), 'r')
            for line in fh:
                line = line.strip().split('\t')
                if 'Pos' not in line: #skip header
                    PSSMs[filename].append({'A' : float(line[1]), 'C' : float(line[2]), 'G' : float(line[3]), 'U' : float(line[4])})

            fh.close()

    return PSSMs

def getkmerscore(kmer, PSSM):
    kmerlen = len(kmer)
    startpos = 0
    kmerscores = []
    PSSMlength = len(PSSM.values()[0]) #number of positions in PSSM
    while startpos <= (PSSMlength - kmerlen):
        ntscores = [] #reset list
        for freqs in PSSM.values(): #slide kmer over 1 by 1, fitting the kmer in the PSSM as many times as you can
            for i in range(0, kmerlen):
                ntscores.append(freqs[i + startpos][kmer[i]])
            kmerscore = reduce(operator.mul, ntscores, 1) #kmerscore is all ntscores multiplied together
            kmerscores.append(kmerscore)
            
            startpos +=1
    
    return max(kmerscores) #return the highest score

def getbestproteinmatch(kmer, PSSMs, masterfile):
    #Return top 1% of protein matches to a kmer as a list of protein names
    
    #Masterfile is SupplementaryData1_RNAcompete_master_file.txt
    #Contains PSSM to protein identities
    topproteins = []
    motiftoprotein = {}
    proteinmatches = {} #{protein : motifscore}

    masterfh = open(masterfile, 'r')
    for line in masterfh:
        line = line.strip().split('\t')
        motifid = line[0]
        proteinid = line[2]
        species = line[3]
        motiftoprotein[motifid] = [proteinid, species]
    masterfh.close()

    kmerscores = {} # {motifID : kmerscore}
    for motif in PSSMs:
        PSSM = {motif : PSSMs[motif]}
        kmerscore = getkmerscore(kmer, PSSM)
        kmerscores[motif] = kmerscore

    #Change motifids to proteins, and only consider human proteins
    for PSSM in kmerscores:
        motifid = PSSM.split('_')[0]
        proteinid = motiftoprotein[motifid][0]
        species = motiftoprotein[motifid][1]
        if species == 'Homo_sapiens': #only consider human proteins
            if proteinmatches.has_key(proteinid) == False:
                proteinmatches[proteinid] = kmerscores[PSSM]

            #Take the maximum score from the two biological replicates (setA and setB)
            elif proteinmatches.has_key(proteinid) == True:
                currentscore = proteinmatches[proteinid]
                newscore = kmerscores[PSSM]
                proteinmatches[proteinid] = max(currentscore, newscore)

    numberofmatches = float(len(proteinmatches))
    top2percent = int(round(numberofmatches / 100) * 2)
    topmatches = sorted(proteinmatches.iteritems(), key = operator.itemgetter(1), reverse = True)[0:top2percent]
    for protein in topmatches:
        topproteins.append(protein[0])

    return topproteins

   

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pssmdirectory', type = str, help = 'Directory containing PSSMs', required = True)
    parser.add_argument('--kmer', type = str, help = 'kmer of interest', required = True)
    parser.add_argument('--masterfile', type = str, help = 'Master file containing motif ID to protein ID conversions.'
                        + 'Usually SupplementaryData1_RNAcompete_master_file.txt', required = True)
    args = parser.parse_args()

    PSSMs= getPSSMs(args.pssmdirectory)
    topproteins = getbestproteinmatch(args.kmer, PSSMs, args.masterfile)
    print topproteins

