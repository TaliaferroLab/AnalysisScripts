#Calculates kmer occurence across a metagene of fasta sequences.  Can handle multiple kmers at once and combine
#their occurences or keep them separate.  Metagene is created as ~50 (exact number depends on length of kmer)
#bins spaced across each sequence in a fasta file.
#Bin width can be changed as indicated in script.
#Counts are normalized so that each bin is counted fairly (some bins get more total kmers than others due to rounding)
#Output is therefore fraction of all kmers at a particular position that are kmer of interest.

#Usage: python kmerpositions.py --help

import sys
import os
import argparse
from Bio import SeqIO

def kmerpositions(kmers, fasta, numberofbins, combinekmers):
    #kmers is list of kmers of interest
    #kmers must all be same length
    #fasta is fasta file of seqs to search through
    #Returns dictionary of kmer counts in ~ 50 bins 
    #(depends on --numberofbins) across sequence 
    binfactor = 100 / float(numberofbins)
    numberofseqs = 0 
    k = len(kmers[0])
    kmerpositionsdict = {} #{binnumber: [{kmer1: counts}, {kmer2: counts}]}
    combinedkmersdict = {} #like kmerpositions dict but counts for all kmers are combined into one value
    normalizationdict = {} #{binnumber: number of positions that feed into that bin} Needed to normalize counts in each bin.
    normalizeddict = {} #like kmerpositionsdict but normalized by the number of positions that feed into every bin

    #intitalize dictionary with 0s for every kmer in every bin
    for record in SeqIO.parse(fasta, 'fasta'):
        numberofseqs +=1
        seqlength = len(record.seq)
        for i in range(seqlength - k + 1):
            position = i + 1
            binnumber = round((float(position / float(seqlength))/binfactor), 2) * binfactor
            kmerpositionsdict[binnumber] = []
            for kmer in kmers:
                kmerdict = {}
                kmerdict[kmer] = 0
                kmerpositionsdict[binnumber].append(kmerdict)
            #Populate normalization dict
            if normalizationdict.has_key(binnumber) == False:
                normalizationdict[binnumber] = 1
            elif normalizationdict.has_key(binnumber):
                normalizationdict[binnumber] += 1

    #Go through sequences, count kmers in each bin
    for record in SeqIO.parse(fasta, 'fasta'):
        seq = str(record.seq.upper().transcribe())
        seqlength = len(record.seq)

        for i in range(seqlength - k + 1):
            kmer = seq[i:i+k]
            position = i + 1
            binnumber = round((float(position / float(seqlength))/binfactor), 2) * binfactor
            if kmer in kmers:
                for item in kmerpositionsdict[binnumber]:
                    if item.keys()[0] == kmer:
                        item[kmer] += 1

    #Normalize the number of counts in each bin by the number of positions that feed into that bin
    for binnumber in kmerpositionsdict:
        normalizeddict[binnumber] = []
        for item in kmerpositionsdict[binnumber]:
            kmer = item.keys()[0]
            counts = item.values()[0]
            normalizedcounts = counts / float(normalizationdict[binnumber])
            kmerdict = {}
            kmerdict[kmer] = normalizedcounts
            normalizeddict[binnumber].append(kmerdict)
            
    #If not combining kmers, you're done.
    if combinekmers == False:
        return normalizeddict
    #If you are combining kmers...
    elif combinekmers == True:
        #blank normalizeddict...gonna remake it
        normalizeddict = {}
        #Sum hits for all kmers in each binnumber
        for binnumber in kmerpositionsdict:
            kmerhits = 0
            for item in kmerpositionsdict[binnumber]:
                kmerhits += item.values()[0]
            combinedkmersdict[binnumber] = kmerhits
        #Normalize counts in each bin by the number fo positions that feed into that bin
        for binnumber in combinedkmersdict:
            normalizedcounts = combinedkmersdict[binnumber] / float(normalizationdict[binnumber])
            normalizeddict[binnumber] = normalizedcounts
        return normalizeddict
    

    
                
                


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type = str, help = 'fasta file of sequences to look through')
    parser.add_argument('--kmers', type = str, help = 'Comma separated (no spaces) list of kmers of interest. Kmers must be of same length. Alternatively, may be single kmer.')
    parser.add_argument('--numberofbins', type = int, help = 'Number of bins in metagene.')
    parser.add_argument('--combinekmers', type = str, help = 'Either <yes> or <no>. Combines all occurences of kmers into one output. Useful when you have mutiple similar kmers for one protein. If only providing one kmer, set to <no>') 
    parser.add_argument('--output', type = str, help = 'Output file.')
    parser.add_argument('--LRbin', type = str, help = 'LR bin number.')
    args = parser.parse_args()

    kmers = args.kmers.upper().replace('T','U').split(',')
    kmerlength = len(kmers[0])
    for kmer in kmers:
        if len(kmer) != kmerlength:
            sys.stderr.write('ERROR: kmers must all be of same length!\n')
            sys.exit()

    if args.combinekmers == 'yes' or args.combinekmers == 'Yes':
        combinekmers = True
        normalizeddict = kmerpositions(kmers, args.fasta, args.numberofbins, combinekmers)
        outfh = open(args.output, 'w')
        outfh.write('bin' + '\t' + 'normalized_kmer_count' + '\t' + 'LRBin' + '\n')
        for binnumber in sorted(normalizeddict):
            outfh.write(str(binnumber) + '\t' + str(normalizeddict[binnumber]) + '\t' + str(args.LRbin) + '\n')
        outfh.close()
    elif args.combinekmers == 'no' or args.combinekmers == 'No':
        combinekmers = False
        normalizeddict = kmerpositions(kmers, args.fasta, args.numberofbins, combinekmers)
        outfh = open(args.output, 'w')
        outfh.write('LRbin' + '\t' + 'bin' + '\t')
        for kmer in kmers:
            outfh.write('normalized_'+kmer+'_count' + '\t')
        outfh.write('\n')
        for binnumber in sorted(normalizeddict):
            outfh.write(str(args.LRbin) + '\t' + str(binnumber))
            for kmer in normalizeddict[binnumber]:
                outfh.write('\t' + str(kmer.values()[0]))
            outfh.write('\n')
        outfh.close()
        
    else:
        sys.stderr.write('ERROR: invalid argument for --combinekmers\n')
        sys.exit()
