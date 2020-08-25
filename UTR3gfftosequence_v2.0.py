#Takes a gff file of 3' UTRs (can have multiple exons) and a directory of fasta files of each chromosome and retrieves sequence for each entry.  Introns are removed and exons are joined together.  This script also expects the stop codon to be already removed.  Also removes last 50 nt of UTR.

#Usage: python 3UTRgfftosequence_v2.0.py <gfffile> <directory_containing_chromosome_sequences> <output.fasta>

#Returns dictionary where key is sequence name and value is sequence.

import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import gffutils

def gfftosequence(gff, sequencedir):

    seqdirectory = os.path.abspath(sequencedir)
    seqfiles = [os.path.join(seqdirectory, file) for file in os.listdir(seqdirectory)]

    #Make index of fasta files
    print 'Indexing sequences...'
    seqdb = SeqIO.index_db('seqdb.idx', seqfiles, 'fasta')
    print '{0} sequences indexed'.format(len(seqdb))

    #Make gff database
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)

    seqs = {} #dictionary where key is ID and value is sequence
    counter = 0

    UTRs = db.features_of_type('mRNA')

    for UTR in UTRs:
        counter +=1
        chrm = UTR.chrom
        strand = UTR.strand
        ID = UTR.attributes['ID']
        exoncoords = []
        UTRsequence = ''
        for exon in db.children(UTR, featuretype = 'exon'):
            exoncoords.append([exon.start, exon.stop])

        number_of_exons = len(exoncoords)

        
        if strand == '+':
            for idx, exonstartstop in enumerate(exoncoords):
                if idx + 1 < number_of_exons:
                    UTRsequence += seqdb[chrm].seq[exonstartstop[0]-1:exonstartstop[1]].upper()
                elif idx + 1 == number_of_exons: #if at the last exon
                    #To check for start codon, take three more nucleotides (change [1]-1 to [1]+2)
                    UTRsequence += seqdb[chrm].seq[exonstartstop[0]-1:exonstartstop[1]-1].upper()

        elif strand == '-':
            for idx, exonstartstop in enumerate(reversed(exoncoords)): #reverse exon order since this is - strand
                if idx + 1 < number_of_exons:
                    UTRsequence += seqdb[chrm].seq[exonstartstop[0]-1:exonstartstop[1]].upper().reverse_complement()
                elif idx + 1 == number_of_exons: #if at the last exon
                    #To check for start codon, take three more nucleotides (change [0] to [0]-3)
                    UTRsequence += seqdb[chrm].seq[exonstartstop[0]:exonstartstop[1]].upper().reverse_complement()
                    
        if len(UTRsequence) > 50:
            seqs[ID] = UTRsequence[:-50] #Remove last 50 nt of 3' UTR
                
        if counter <= 50 and counter % 10 == 0:
            print 'Retrieving sequence %i' % (counter)
        elif counter > 50 and counter % 50 == 0:
            print 'Retrieving sequence %i' % (counter)

    print 'Retrieved {0} sequences.'.format(len(seqs))

    os.remove(db_fn)
    
    os.remove('seqdb.idx')
    
    return seqs
    
if __name__ == '__main__':
    seqdictionary = gfftosequence(sys.argv[1], sys.argv[2])
    outfh = open(sys.argv[3], 'w')
    for ID in seqdictionary:
        outfh.write('>' + ID + '\n' + str(seqdictionary[ID]) + '\n') #print in fasta format

    outfh.close()


    
