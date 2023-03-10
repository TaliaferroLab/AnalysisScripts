#Given a sequence, chop it into pieces of length n every x bp.
#Then, fold those sequences and look for the bp probability of specific kmers within those sequences

import subprocess
from Bio import SeqIO
import argparse
import re
import os
import pandas as pd
import sys


def breaksequence(seq, chunksize, chunkstep):
    #Given a sequence, break it into chunks of size chunksize with one chunk every chunkstep
    chunkseqs = []
    for i in range(0, len(seq), chunkstep):
        if i + chunksize <= len(seq):
            chunkseq = seq[i : i + chunksize]
            chunkseqs.append(chunkseq)
            i += chunkstep
        else:
            break

    #If the last chunk does not lie flush with the end of the sequence, add one more that does so that the end is always covered
    if (i - chunkstep) + chunksize != len(seq):
        chunkseq = seq[chunksize * -1 :]
        chunkseqs.append(chunkseq) 

    return chunkseqs

def getkmerpositions(seq, kmer):
    #Get kmer positions within a sequence
    mp = [m.start() for m in re.finditer('(?={0})'.format(kmer), seq)]
    motifcount = len(mp)
    kmerlength = len(kmer)
    # 0-indexed list of all motif positions (any position covered by kmer)
    motifpositions = []

    for motifpos in mp:
        motifoccurence = []
        for i in range(kmerlength):
            motifoccurence.append(motifpos + i)
        motifpositions.append(motifoccurence)

    return motifpositions, motifcount

def getsummedbpprobs(seq, kmer):
    #Given a sequence, make bp probabilities and return the summed bp probability for a base (the sum of its pairing to every other base)
    probs = {}  # {base : bpprob}
    motifprobs = {} #{motifposition1 : [probs], motifposition2 : [probs]}

    motifpositions, motifcount = getkmerpositions(seq, kmer)
    

    #Populate dictionary
    for i in range(len(seq)):
        probs[i] = 0
    for i in range(len(kmer)):
        motifprobs[i] = []

    command = 'RNAfold -p'
    job = subprocess.Popen(command, shell = True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    job.stdin.write(seq.encode('utf-8'))
    output = job.communicate()
    #structure = str(output[0].split('\n')[1].split(' ')[0]).encode('utf-8')
    #print structure
    bpfh = open('dot.ps', 'r')
    #lines containing bp probs have 'ubox' in the 4th field
    for line in bpfh:
        line = line.strip().split(' ')
        if len(line) != 4:
            continue
        if line[3] == 'ubox':
            leftbase = int(line[0]) - 1  # make these 0 based
            rightbase = int(line[1]) - 1
            #3rd field is square root of bpprob
            bpprob = float(line[2])**2
            probs[leftbase] += bpprob
            probs[rightbase] += bpprob

    bpfh.close()
    os.remove('dot.ps')
    os.remove('rna.ps')

    #Get the bp prob over the motif
    for motifoccurence in motifpositions:
        for idx, i in enumerate(motifoccurence):
            motifprobs[idx].append(round(probs[i], 3))

    return motifprobs, motifcount

def iteratefasta(fasta, output, chunksize, chunkstep, kmer):
    #Given a fasta of sequences of interest, for each sequence,
    #(1), break it into chunks with breaksequence()
    #(2), for each chunk, fold it and record bp probs over motif (if motif exists in chunk)
    #(3), record motif bp probs in table:
    #one row per motif occurence (motifs can end up in multiple rows if they show up in multiple chunks of a given fasta seq)
    #one column for sequence name (from fasta), one column for each motif position

    outdict = {} #{seqname = [list of seqs interrogated (can be repeats from chunks of same seq)], motifposition1 = [bpprobs], motifposition2[bpprobs], ...}
    kmer = kmer.replace('T', 'U').upper()

    #Initialize dictionary
    outdict['seqname'] = []
    for i in range(len(kmer)):
        outdict[i] = []
    
    seqcounter = 0
    for record in SeqIO.parse(fasta, 'fasta'):
        seqcounter +=1
        if seqcounter % 1000 == 0:
            print('Folding sequence {0}...'.format(seqcounter))
        recordid = str(record.id)
        seq = str(record.seq).replace('T', 'U').upper()
        #REMOVING ADAPTERS FROM OLIGOS....TAKE THIS OUT IF NOT NEEDED
        seq = seq[20:-20]
        chunkseqs = breaksequence(seq, int(chunksize), int(chunkstep))
        for chunkseq in chunkseqs:
            #no use folding this if there isn't a kmer in there
            if kmer in chunkseq:
                motifprobs, motifcount = getsummedbpprobs(chunkseq, kmer)
                for i in range(motifcount):
                    outdict['seqname'].append(recordid)
                for i in range(len(kmer)):
                    outdict[i].extend(motifprobs[i])

    #Turn into df for output
    df = pd.DataFrame.from_dict(outdict)
    df.to_csv(output, sep = '\t', index = False)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required = True, type = str, help = 'Fasta file of sequences to fold.')
    parser.add_argument('--kmer', required = True, type = str, help = 'kmer to consider for bp probabilities.')
    parser.add_argument('--chunksize', required = True, type = int, help = 'Size of sequence chunks for folding to be make from each fasta entry.')
    parser.add_argument('--chunkstep', required = True, type = int, help = 'Step size between sequence chunks.')
    parser.add_argument('--output', required = True, help = 'Output file.')
    args = parser.parse_args()

    iteratefasta(args.fasta, args.output, args.chunksize, args.chunkstep, args.kmer)