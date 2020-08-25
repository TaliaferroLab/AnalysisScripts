#Take the output files from shape_normalizereactivites.py.  These include a directory of fasta files for each sequence,
#a directory of SHAPE files for each sequence, and the normalized reactivities.out.  Run RNAstructure software for each
#sequence to get a basepair probability value for each nt.  Append that value to the normalized reactivities.out file
#for each nt.


import argparse
import subprocess
import os
import sys


#The output from ProbabilityPlot is the probability of all possible basepairs. I just want to know
#what's the probability of a base to be paired to any other base.  So, for each nucleotide, sum all of
#the basepair probabilities for that nucleotide.
def summarizebpprobmatrix(bpprobmatrix):
    infh = open(bpprobmatrix, 'r')
    bpprobs = {} # {nt : cumulative_bpprob}
    for line in infh:
        line = line.strip().split('\t')
        if len(line) != 3 or line[0] == 'i': #skip header lines
            continue
        nt1 = int(line[0])
        nt2 = int(line[1])
        prob = 10**(float(line[2]) * -1)
        if nt1 not in bpprobs:
            bpprobs[nt1] = prob
        elif nt1 in bpprobs:
            newprob = bpprobs[nt1] + prob
            bpprobs[nt1] = newprob
        if nt2 not in bpprobs:
            bpprobs[nt2] = prob
        elif nt2 in bpprobs:
            newprob = bpprobs[nt2] + prob
            bpprobs[nt2] = newprob

    infh.close()

    return bpprobs
        


#Given a directory of fasta files and SHAPE outputs (from shape_normalizereactivities.py), use RNAstructure
#to calculate partition functions and then basepair probabilities for every nt. The basepair probabilities output
#by ProbabilityPlot are matricies of all possible basepairs and should be summarized using summarizebpprobmatrix.
def getbpprobs(fastadir, shapedir):
    fastadir = os.path.abspath(fastadir)
    shapedir = os.path.abspath(shapedir)
    allbpprobs = {} # {seqname : {nt : bpprob}}
    for fasta in os.listdir(fastadir):
        seq = os.path.basename(fasta)
        fastafile = os.path.join(fastadir, seq)
        shapefile = os.path.join(shapedir, seq[:-5] + 'SHAPE')
        print 'Calculating partition function for {0}...'.format(seq)
        if os.path.isfile(shapefile):
            subprocess.call(['partition', fastafile, '--SHAPE', shapefile, 'output.pfs'])
        elif os.path.isfile(shapefile) == False:
            #if there's not a shape file, that's becuase it didn't pass the --minstops filter when the 
            #SHAPE file was created with shape_normalize_reactivities.py
            print 'No SHAPE file found for {0}!!!'.format(seq)
            continue
        print 'Calculating basepair probability matrix for {0}...'.format(seq)
        subprocess.call(['ProbabilityPlot', '--text', 'output.pfs', 'bpprobmatrix.txt'])
        bpprobs = summarizebpprobmatrix('bpprobmatrix.txt')
        allbpprobs[seq[:-6]] = bpprobs
        os.remove('output.pfs')
        os.remove('bpprobmatrix.txt')

    return allbpprobs

def addbpprobstoreactivities(bpprobs, reactivities, outfile):
    outfh = open(outfile, 'w')
    infh = open(reactivities, 'r')
    for line in infh:
        line = line.strip().split('\t')
        if line[0] == 'sequence': #if this is the header line
            outfh.write(('\t').join(line) + '\t' + 'bp_prob' + '\n')
        else:
            seqname = line[0]
            nt = int(line[2])
            if seqname in bpprobs: 
                #only write those lines for which RNAstructure calculated bp probs
                #this corresponds to those RNAs that passed the --minstops filter when the SHAPE file was
                #created with shape_normalizereactivities.py
                if nt == 0:
                    outfh.write(('\t').join(line) + '\t' + '0' + '\n')
                elif nt > 0:
                    bp_prob = str(bpprobs[seqname][nt])
                    outfh.write(('\t').join(line) + '\t' + bp_prob + '\n')

            else:
                continue

    infh.close()
    outfh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastadir', type = str, help = 'Directory of fasta sequences ouput from shape_normalizereactivities.py.')
    parser.add_argument('--shapedir', type = str, help = 'Directory of shape files ouput from shape_normalizereactivities.py.')
    parser.add_argument('--reactivities', type = str, help = 'Normalized reactivities output from shape_normalizereactivities.py.')
    parser.add_argument('--outfile', type = str, help = 'Output file. Will look exactly like reactivities file but with bp probabilities field added on.')
    args = parser.parse_args()

    bpprobs = getbpprobs(args.fastadir, args.shapedir)
    addbpprobstoreactivities(bpprobs, args.reactivities, args.outfile)
