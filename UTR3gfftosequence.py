#Takes a gff file and a directory of fasta files for each chromosome and retrieves sequence for each entry.

#Usage: python gfftosequence.py <gfffile> <directory_containing_chromosome_sequences> <output.fasta>

#Returns dictionary where key is sequence name and value is sequence.  If used standalone, prints output in fasta format.

import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np

def gfftosequence(gff, sequencedir):

    seqdirectory = os.path.abspath(sequencedir)
    seqfiles = [os.path.join(seqdirectory, file) for file in os.listdir(seqdirectory)]

    seqdb = SeqIO.index_db('seqdb.idx', seqfiles, 'fasta') #make index of fasta files
    print '%i sequences indexed' % len(seqdb)
    
    gfffh = open(gff, 'r')
    gffentries = []
    seqs = {} #dictionary where key is ID and value is sequence
    GCs = []
    counter = 0

    for line in gfffh:
        line = line.strip().replace(';','\t').split('\t')
        gffentries.append(line)
    gfffh.close()

    for gffentry in gffentries:
        chrm = gffentry[0]
        start = int(gffentry[3])
        stop = int(gffentry[4])
        strand = gffentry[6]
        ID = gffentry[8][3:] #expecting this field to be ID=...

        if strand == '+':
            sequence = seqdb[chrm].seq[start+2:stop+1].upper() #can slice differently here...this removes stop codons
        elif strand == '-':
            sequence = seqdb[chrm].seq[start-1:stop-3].reverse_complement().upper()
        seqs[ID] = sequence

        GCs.append(float(GC(sequence)))

        counter +=1
        if counter <= 50 and counter % 10 == 0:
            print 'Retrieving sequence %i of %i' % (counter,len(gffentries))
        elif counter > 50 and counter % 50 == 0:
            print 'Retrieving sequence %i of %i' % (counter,len(gffentries))

    print 'The average GC content of these sequences is %f' % np.mean(GCs)
    os.remove('seqdb.idx')

    return seqs 
       

if __name__ == '__main__':
    seqdictionary = gfftosequence(sys.argv[1], sys.argv[2])
    outfh = open(sys.argv[3], 'w')

    for ID in seqdictionary:
        outfh.write('>' + ID + '\n' + str(seqdictionary[ID]) + '\n') #print in fasta format

    outfh.close()
