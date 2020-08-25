#Phastcons.bed is 0-based.  UTR gff is 1-based.
#Necessary modules: biopython, pysam
#Usage: python motifconservations_phastcons.py <UTR_gff> <genome.fasta> <k> <sortedphastconsbed.gz> <kmer_list (output of kmerenrichment.py> <simple/complex> <number_of_control_iterations>
#Must also have sortedphastconsbed.gz.tbi in same directory as sortedphastconsbed.gz

#Not sure if this is a good idea.  Probably better to use motifconservation_v2.0.py

import sys
import os
import random
import gffutils
import tabix
import pysam
from Bio import SeqIO
import numpy as np

def getSequences(gff, sequence_file):
    sys.stderr.write('Indexing genome sequences...\n')
    seq_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, 'fasta'))
    sys.stderr.write('Indexed {0} sequences.\n'.format(len(seq_dict)))

    #Make gff database
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)
    db = gffutils.FeatureDB(db_fn)
    seqs = {} #fasta format in dictionary {seqID:sequence}
    UTRs = db.features_of_type('3\'UTR')

    sys.stderr.write('Retrieving sequences from gff...\n')
    for UTR in UTRs:
        chrm = str(UTR.chrom)
        strand = str(UTR.strand)
        start = int(UTR.start)
        stop = int(UTR.stop)
        ID = str(UTR.id)
        if strand == '+':
            seqID = ID + ';' + chrm + ';' + str(start+3) + ';' + str(stop-50) + ';' + strand #gonna take off stop codon and last 50 nt
        elif strand == '-':
            seqID = ID + ';' + chrm + ';' + str(start+50) + ';' + str(stop-3) + ';' + strand #gonna take off stop codon and last 50 nt
        UTRsequence = ''

        if strand == '+':
            UTRsequence += seq_dict[chrm].seq[start+2:stop].upper() #remove stop codons
        elif strand == '-':
            UTRsequence += seq_dict[chrm].seq[start-1:stop-3].upper().reverse_complement() #remove stop codons

        if len(UTRsequence) > 50:
            seqs[seqID] = str(UTRsequence[:-50])

    os.remove(db_fn)
    sys.stderr.write('Retrieved {0} UTR sequences over 50 nt.\n'.format(len(seqs)))
    return seqs 

def getKmerPhastcons(fastadict, k, phastconsbed):
    k = int(k)
    phastconsdict = {} #{kmer:mean_phastcons_score_of_bases_in_kmer}
    phastconsaveragedict = {} #{kmer:mean of scores in phastconsdict}
    UTRcounter = 0
    phastconstabix = pysam.Tabixfile(phastconsbed)

    for UTR in fastadict:
        UTRcounter +=1
        if UTRcounter % 50 == 0:
            sys.stderr.write('Determining motif conservation in UTR {0} of {1}...\n'.format(UTRcounter, len(fastadict)))
        UTRsequence = fastadict[UTR]
        UTR = UTR.replace(';', '\t').split('\t')
        ID = UTR[0]
        chrm = UTR[1]
        start = int(UTR[2])
        stop = int(UTR[3])
        strand = UTR[4]
        for i in range(len(UTRsequence) - k + 1):
            if strand == '+':
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = start + i
                mousekmerstop = start + i + k - 1
            elif strand == '-':
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = stop - i - k + 1
                mousekmerstop = stop - i
            kmerscores = [] #list of phastcons scores for every bp of this kmer
            #Remember...input gff coords are 1-based and the phastconsbed is 0-based
            for bed in phastconstabix.fetch(str(chrm), mousekmerstart - 1, mousekmerstop, parser = pysam.asBed()):
                kmerscores.append(float(bed.name))
            if len(kmerscores) == k: #if every base in the kmer had a score
                kmeraveragescore = (sum(kmerscores) / float(len(kmerscores)))
                if phastconsdict.has_key(mousekmer) == False:
                    phastconsdict[mousekmer] = [kmeraveragescore] #if kmer not in dictionary, initialize entry
                elif phastconsdict.has_key(mousekmer):
                    phastconsdict[mousekmer].append(kmeraveragescore)
            elif len(kmerscores) != k: #if not every base in the kmer had a score
                continue

    for kmer in phastconsdict:
        phastconsaveragedict[kmer] = np.mean(phastconsdict[kmer])
    return phastconsaveragedict

#Use if kmer is not complex.  Returns random kmers with same GC and CpG contents.
def randomkmer(kmer, iterations):
    iterations = int(iterations)
    sys.stderr.write('Creating {0} random kmers with matched GC content and CpG counts for kmer {1}...\n'.format(iterations, kmer))
    controlsequences = []
    kmer = kmer.upper()
    CpGcount = kmer.count('CG')
    GCcontent = kmer.count('C') + kmer.count('G')
    kmersize = len(kmer)
    while len(controlsequences) < iterations:
        randomseq = ''.join(random.choice('ACGT') for x in range(kmersize))
        randomCpGcount = randomseq.count('CG')
        randomGCcontent = randomseq.count('C') + randomseq.count('G')
        if randomCpGcount == CpGcount and randomGCcontent == GCcontent and randomseq != kmer:
            controlsequences.append(randomseq)
            
    sys.stderr.write('Done!\n')
    return controlsequences 

#Use if kmer is complex. Returns shuffled versions of that kmer with matching CpG contents.
def shufflekmer(kmer, iterations):
    iterations = int(iterations)
    sys.stderr.write('Shuffling kmer {0} times and matching CpG count to create control kmers... \n'.format(iterations))
    controlsequences = []
    kmer = kmer.upper()
    CpGcount = kmer.count('CG')
    kmer = list(kmer)
    while len(controlsequences) < iterations:
        random.shuffle(kmer)
        shuffled = ''.join(kmer)
        if shuffled.count('CG') == CpGcount and shuffled != kmer:
            controlsequences.append(shuffled)

    sys.stderr.write('Done!\n')
    return controlsequences

def summarizePhastconsAverageDict(phastconsaveragedict, kmer_of_interest, kmerlist):
    controlkmer_phastconsscores = []
    if phastconsaveragedict.has_key(kmer_of_interest) == False:
        kmer_phastconsscore = 'NA'
        controlkmer_mean = 'NA'
        controlkmer_std = 'NA'
        zscore = 'NA'
    elif phastconsaveragedict.has_key(kmer_of_interest):
        for kmer in phastconsaveragedict:
            if kmer != kmer_of_interest and kmer in kmerlist:
                controlkmer_phastconsscores.append(phastconsaveragedict[kmer])
            elif kmer == kmer_of_interest and kmer in kmerlist:
                kmer_phastconsscore = phastconsaveragedict[kmer]

        controlkmer_mean = np.mean(controlkmer_phastconsscores)
        controlkmer_std = np.std(controlkmer_phastconsscores)
        if controlkmer_std != 0:
            zscore = ((kmer_phastconsscore - controlkmer_mean) / float(controlkmer_std))
        elif controlkmer_std == 0:
            zscore = 'NA because control std was 0'

    kmerstats = [str(kmer_phastconsscore), str(controlkmer_mean), str(controlkmer_std), str(zscore)]

    return kmerstats
             
                
            
if __name__ == '__main__':
    fastadict = getSequences(sys.argv[1], sys.argv[2])
    phastconsaveragedict = getKmerPhastcons(fastadict, sys.argv[3], sys.argv[4])
    kmerlistfile = open(sys.argv[5], 'r')
    kmer_complexity = sys.argv[6]
    control_iterations = sys.argv[7]
    fh = open('kmerconservations_phastcons.txt', 'w')
    fh.close()
    for line in kmerlistfile:
        line = line.strip().split('\t')
        if line[0] != 'kmer':
            kmer_of_interest = line[0].upper().replace('U','T')
            if kmer_of_interest == 'CGCGCG': #impossible to get control kmers for this
                continue
            kmerenrichment = str(line[3])
            pvalue = str(line[4])
            bhpvalue = str(line[5])
            if kmer_complexity == 'simple':
                kmerlist = randomkmer(kmer_of_interest, control_iterations) + [kmer_of_interest]
            elif kmer_complexity == 'complex':
                kmerlist = shufflekmer(kmer_of_interest, control_iterations) + [kmer_of_interest]

            #Summarize with summarizePhastconsAverageDict
            kmerstats = summarizePhastconsAverageDict(phastconsaveragedict, kmer_of_interest, kmerlist)
            outfh = open('kmerconservations_phastcons.txt', 'a')
            outfh.write(kmer_of_interest.replace('T','U') + '\t' + line[1] + '\t' + line[2] + '\t' + kmerenrichment + '\t' + pvalue + '\t' + bhpvalue + '\t' + kmerstats[0] + '\t' + kmerstats[1] + '\t' + kmerstats[2] + '\t' + kmerstats[3] + '\n')
            outfh.close()
    kmerlistfile.close()
    sys.stderr.write('All done!\n')

