from Bio import SeqIO
from itertools import izip
import itertools
import argparse
from Bio.Alphabet import IUPAC

#First get all sequences that have a particular barcode
#Barcode is first four nts on reverse read
#Allows one mismatch in barcodes
def getBarcodedSeqs(barcode1, barcode2, fastqread1, fastqread2, forwardout1, reverseout1, forwardout2, reverseout2):
    currentseq = 0
    currentrevseq = 0
    counter = 0
    forwardreads = []
    reversereads = []
    fastq_iter1 = SeqIO.parse(fastqread1, 'fastq')
    fastq_iter2 = SeqIO.parse(fastqread2, 'fastq')
    for currentseq, (rec1, rec2) in enumerate(itertools.izip(fastq_iter1, fastq_iter2)):
        counter +=1
        if counter % 1000000 == 0:
            print 'Processing read {0}'.format(counter)

        #For shape v1.0 with degenerate barcode on forward read
        seqbarcode = str(rec1.seq[0:4]).replace('T','Y').replace('C','Y').replace('A','R').replace('G','R')
        #For shape v2.0 with barcode on reverse read
        #seqbarcode = rec2.seq[0:4]
        if sum(ii == jj for ii, jj in izip(seqbarcode, barcode1)) >= 3:
            with open(forwardout1, 'a') as forwardreads:
                SeqIO.write(rec1, forwardreads, 'fastq')

            #Get rid of any bases (up to 3 of them) upstream of an obvious instance where the RNA is full length
            #This only works if the 5' end of every RNA is the same (here it begins with GGGCCTTGACACC), which was
            #true for the original ancient exon pool that had Bind-n-seq adapters on both ends. 
            if str(rec2.seq[1:14]) == 'GGGCCTTGACACC':
                with open(reverseout1, 'a') as reversereads:
                    SeqIO.write(rec2[1:], reversereads, 'fastq')

            elif str(rec2.seq[2:15]) == 'GGGCCTTGACACC':
                with open(reverseout1, 'a') as reversereads:
                    SeqIO.write(rec2[2:], reversereads, 'fastq')

            elif str(rec2.seq[3:16]) == 'GGGCCTTGACACC':
                with open(reverseout1, 'a') as reversereads:
                    SeqIO.write(rec2[3:], reversereads, 'fastq')

            else:
                with open(reverseout1, 'a') as reversereads:
                    SeqIO.write(rec2, reversereads, 'fastq')

        elif sum(ii == jj for ii, jj in izip(seqbarcode, barcode2)) >=3:
            with open(forwardout2, 'a') as forwardreads:
                SeqIO.write(rec1, forwardreads, 'fastq')

            #Get rid of any bases (up to 3 of them) upstream of an obvious instance where the RNA is full length
            if str(rec2.seq[1:14]) == 'GGGCCTTGACACC':
                with open(reverseout2, 'a') as reversereads:
                    SeqIO.write(rec2[1:], reversereads, 'fastq')

            elif str(rec2.seq[2:15]) == 'GGGCCTTGACACC':
                with open(reverseout2, 'a') as reversereads:
                    SeqIO.write(rec2[2:], reversereads, 'fastq')

            elif str(rec2.seq[3:16]) == 'GGGCCTTGACACC':
                with open(reverseout2, 'a') as reversereads:
                    SeqIO.write(rec2[3:], reversereads, 'fastq')

            else:
                with open(reverseout2, 'a') as reversereads:
                    SeqIO.write(rec2, reversereads, 'fastq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode1', type = str, help = 'First barcode to look for. Should be first 4 nts reverse read.', required = True)
    parser.add_argument('--barcode2', type = str, help = 'Second barcode to look for. Should be first 4 nts reverse read.', required = True)
    parser.add_argument('--forwardfastq', type = str, help = 'Fastq file of first seqeuncing read.', required = True)
    parser.add_argument('--reversefastq', type = str, help = 'Fastq file of second sequencing read.', required = True)
    parser.add_argument('--forwardout1', type = str, help = 'Outfile for forward reads with barcode 1.', required = True)
    parser.add_argument('--reverseout1', type = str, help = 'Outfile for reverse reads with barcode 1.', required = True)
    parser.add_argument('--forwardout2', type = str, help = 'Outfile for forward reads with barcode 2.', required = True)
    parser.add_argument('--reverseout2', type = str, help = 'Outfile for reverse reads with barcode 2.', required = True)
    args = parser.parse_args()
    forwardout1 = open(args.forwardout1, 'w')
    reverseout1 = open(args.reverseout1, 'w')
    forwardout2 = open(args.forwardout2, 'w')
    reverseout2 = open(args.reverseout2, 'w')
    

    getBarcodedSeqs(args.barcode1, args.barcode2, args.forwardfastq, args.reversefastq, args.forwardout1, args.reverseout1, args.forwardout2, args.reverseout2)

