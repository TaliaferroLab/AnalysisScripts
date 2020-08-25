#Fasta1 = 'test sequences' ; Fasta2 = 'background sequences'
#Usage: python kmerenrichment.py -h
#Returns: <kmer> <fastafile1count> <fastafile2count> <enrichment> <pvalue> <bh_adjusted_pvalue>

import operator
import sys
from Bio import SeqIO
from scipy.stats import fisher_exact
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import argparse
from kmertoprotein import getPSSMs, getbestproteinmatch


def countKmers(fasta1, fasta2, k):
    kmers1 = {}
    kmers2 = {}
    currentkmer = ''
    previouskmer = ''
    currentkmercount = 0
    cont_table = {}
    bh_cont_table = {}
    pvalues = []
    bh_pvalues = []
    output_list = []
    counter = 0
    k = int(k)

    #Count kmers in fasta1
    for seq_record in SeqIO.parse(fasta1, 'fasta'):
        seq = str(seq_record.seq.transcribe())

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            #Count homopolymeric stretches as only one instance of kmer, but allow nonoverlapping instances of
            #kmer in homopolymers
            if ((kmer != currentkmer and kmer != previouskmer) or currentkmercount >= k) and 'N' not in kmer: 
                if kmers1.has_key(kmer):
                    kmers1[kmer] +=1
                else:
                    kmers1[kmer] = 1
                #To avoid over-counting dinucleotide repeats (UCUCUCUCUCUCUCUC, for example), keep track of what the kmer 2 rounds ago was.
                #currentkmer is the kmer from last round, previous kmer is from 2 rounds ago.
                previouskmer = currentkmer
                currentkmer = kmer
                currentkmercount = 1
            else:
                currentkmercount +=1
                
    currentkmer = '' #reset currentkmer
    previouskmer = ''
    currentkmercount = 0
    
    #Count kmers in fasta2    
    for seq_record in SeqIO.parse(fasta2, 'fasta'):
        seq = str(seq_record.seq.transcribe())

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if ((kmer != currentkmer and kmer != previouskmer) or currentkmercount >= k) and 'N' not in kmer:
                if kmers2.has_key(kmer):
                    kmers2[kmer] +=1
                else:
                    kmers2[kmer] = 1
                previouskmer = currentkmer
                currentkmer = kmer
                currentkmercount = 1
            else:
                currentkmercount +=1
                   
    #Make sure all kmers present in one dictionary are present in the other
    for kmer in kmers1:
        if kmer not in kmers2:
            kmers2[kmer] = 0

    for kmer in kmers2:
        if kmer not in kmers1:
            kmers1[kmer] = 0


    #Get total number of kmers in each file
    totalkmers1 = sum(kmers1.values())
    totalkmers2 = sum(kmers2.values())

    #Make contigency table for each kmer
    for kmer in kmers1:
        cont_table[kmer] = [kmers1[kmer], kmers2[kmer], totalkmers1 - kmers1[kmer], totalkmers2 - kmers2[kmer]]

    #Peform fisher's exact test
    for kmer in cont_table:
        fasta1kmercount = cont_table[kmer][0]
        fasta2kmercount = cont_table[kmer][1]
        fasta1otherkmers = cont_table[kmer][2]
        fasta2otherkmers = cont_table[kmer][3]
        oddsratio, pvalue = fisher_exact([[fasta1kmercount, fasta2kmercount], [fasta1otherkmers, fasta2otherkmers]])

        #Calculate enrichment to two decimal places
        if fasta1kmercount > 0 and fasta2kmercount > 0 :
            enrichment = '%.2f' % ((float(fasta1kmercount)/float(totalkmers1)) / (float(fasta2kmercount)/float(totalkmers2)))
            cont_table[kmer].append(enrichment)
        else:
            cont_table[kmer].append('NA')
        cont_table[kmer].append(pvalue)
        counter +=1

        if counter % 100 == 0:
            print 'Analyzing kmer %i of %i' % (counter, len(kmers1))

        
        
    # Sort by pvalue
    sortedcont_table = sorted(cont_table.items(), key = lambda x: x[1][5])

    #Calculate Benjamini-Hochsberg corrected pvalue
    stats = importr('stats')
    for kmer in sortedcont_table:
        pvalues.append(kmer[1][5])

    p_adjust = stats.p_adjust(FloatVector(pvalues), method = 'BH')

    for pvalue in p_adjust:
        bh_pvalues.append(pvalue)

    for kmeridx, kmer in enumerate(sortedcont_table):
        for bh_pvalueidx, bh_pvalue in enumerate(bh_pvalues):
            if kmeridx == bh_pvalueidx:
                bh_cont_table[kmer[0]] = [kmer[1][0], kmer[1][1], kmer[1][2], kmer[1][3], kmer[1][4], kmer[1][5], bh_pvalue]

    sortedbh_cont_table = sorted(bh_cont_table.items(), key = lambda x: x[1][6])

    #Make a list for output
    for kmer in sortedbh_cont_table:
        output_list.append([kmer[0], kmer[1][0], kmer[1][1], kmer[1][4], kmer[1][5], kmer[1][6]])

    
    return output_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences', type = str, help = 'Fasta file of sequences you care about.', required = True)
    parser.add_argument('--background', type = str, help = 'Fasta file of control or background sequences.', required = True)
    parser.add_argument('-k', type = int, help = 'kmer length to search for.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    parser.add_argument('--kmertoprotein', help = 'Use to find the best protein matches for each identified kmer.', action = 'store_true', required = False)
    parser.add_argument('--pssmdirectory', type = str, help = 'Directory containing PSSMs. Required if using kmertoprotein.', required = False)
    parser.add_argument('--masterfile', type = str, help = 'Master file containing motif ID to protein ID conversions.'
                        + 'Usually SupplementaryData1_RNAcompete_master_file.txt. Required if using kmertoprotein.', required = False)
    args = parser.parse_args()
    
    outfh = open(args.outfile, 'w')

    if args.kmertoprotein == False:
        outfh.write('kmer' + '\t' + str(args.sequences) + '_count' + '\t' + str(args.background) + '_count' + '\t' 'enrichment' + '\t' + 'fisher\'s exact p' + '\t' + 'BH corrected p''\n')
        for kmer in countKmers(args.sequences, args.background, args.k):
            outfh.write(('\t').join([str(item) for item in kmer]) + '\n')

    elif args.kmertoprotein == True:
        PSSMs = getPSSMs(args.pssmdirectory)
        outfh.write('kmer' + '\t' + str(args.sequences) + '_count' + '\t' + str(args.background) + '_count' + '\t' 'enrichment' + '\t' + 'fisher\'s exact p' + '\t' + 'BH corrected p' + '\t' + 'Best_protein_matches''\n')
        for kmer in countKmers(args.sequences, args.background, args.k):
            kmerseq = kmer[0]
            topproteinslist = getbestproteinmatch(kmerseq, PSSMs, args.masterfile) #this is a list
            topproteins = (',').join(topproteinslist)
            outfh.write(('\t').join([str(item) for item in kmer]) + '\t' + topproteins + '\n')

    outfh.close()


