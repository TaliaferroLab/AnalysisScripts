from Bio import SeqIO
from itertools import izip
import itertools
import argparse

#Usage: python filterreadlength_paired.py --help

#Requires both a read AND its mate to pass length filter.

def trimSeqs(fastqread1, fastqread2, minlength, forwardout, reverseout):
    minlength = int(minlength)
    currentseq = 0
    currentrevseq = 0
    forwardreads = []
    reversereads = []
    fastq_iter1 = SeqIO.parse(fastqread1, 'fastq')
    fastq_iter2 = SeqIO.parse(fastqread2, 'fastq')
    for currentseq, (rec1, rec2) in enumerate(itertools.izip(fastq_iter1, fastq_iter2)):
        if currentseq % 1000000 == 0:
            print 'Processing read {0}'.format(currentseq)
        if len(rec1) >= minlength and len(rec2) >= minlength:
            with open(forwardout, 'a') as forwardoutfh:
                SeqIO.write(rec1[:minlength], forwardoutfh, 'fastq')
            with open(reverseout, 'a') as reverseoutfh:
                SeqIO.write(rec2[:minlength], reverseoutfh, 'fastq')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--forwardfastq', type = str, help = 'Fastq file of first seqeuncing read.', required = True)
    parser.add_argument('--reversefastq', type = str, help = 'Fastq file of second sequencing read.', required = True)
    parser.add_argument('--forwardout', type = str, help = 'Outfile for forward reads.', required = True)
    parser.add_argument('--reverseout', type = str, help = 'Outfile for reverse reads.', required = True)
    parser.add_argument('--minlength', type = int, help = 'Minimum sequence length to keep. All seqs longer than minlength will be trimmed to minlength.', required = True)
    args = parser.parse_args()
    
    forwardoutfh = open(args.forwardout, 'w')
    reverseoutfh = open(args.reverseout, 'w')
    forwardoutfh.close()
    reverseoutfh.close()
    trimSeqs(args.forwardfastq, args.reversefastq, args.minlength, args.forwardout, args.reverseout)
    
