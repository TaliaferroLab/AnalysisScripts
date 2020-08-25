#Takes a MISO gff of ALE events and a reference gff. Reference must have stop codons explicitly annotated as such.  Searches for stop codons in each ALE.  If none, prints a warning and skips.  If >1, takes the stop codon that is most proximal to the gene.  Returns gff of CDS portion of ALE, the nucleotide sequence of the CDS portion in fasta format, and a translation of that sequence also in fasta format. 
#Only reports translation if there is 1 and only 1 stop codon in translated sequence.

#ALE isoforms must have class 'mRNA' in 3rd field, as is customary for MISO annotations.

#ASSUMES ALE ISOFORMS ARE SINGLE EXONS.  This is true almost all of the time.

#Holds entire genome sequence in memory.

#Standalone usage: python getALECDS.py -h

import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np

def getALEORFCoords(ALEgff, stopcodongff):

    ALEfh = open(ALEgff, 'r')
    stopcodonfh = open(stopcodongff, 'r')
    
    stopcodons = []
    ALEs = []
    
    for line in stopcodonfh:
        line = line.strip()
        line = line.split('\t')
        if line[2] == 'stop_codon':
            stopcodons.append(line)

    for line in ALEfh:
        line = line.strip()
        line = line.split('\t')
        if line[2] == 'mRNA':
            ALEs.append(line)

    ALEfh.close()
    stopcodonfh.close()


    stoplist=[]
    stoppos = []
    ORFgffAll = []

    for ALE in ALEs:
        ALEchrm = ALE[0]
        ALEstart = int(ALE[3])
        ALEstop = int(ALE[4])
        ALEstrand = ALE[6]

        for stopcodon in stopcodons:
            stopcodonchrm = stopcodon[0]
            stopcodonstart = int(stopcodon[3])
            stopcodonstop = int(stopcodon[4])
            stopcodonstrand = stopcodon[6]
            
            #if chromosomes and strands match AND stop codon start/stop lies between ALE start/stop
            if ALEchrm == stopcodonchrm and ALEstrand == stopcodonstrand and ALEstart <= stopcodonstart and ALEstop >= stopcodonstop:
                if ALEstrand == '+': #append beginning of start codon...this will include stop codon itself, useful for checking
                    stoplist.append(stopcodonstart)
                    
                elif ALEstrand == '-': #append end of start codon, which since on the - strand, is the beginning of it
                    stoplist.append(stopcodonstop)
        #Remove duplicates
        print len(list(set(stoplist))) 
        

        if len(stoplist) > 0:
            if ALEstrand == '+': # take most proximal stop codon to give longest UTR and shortest CDS
                stoppos = [ALEchrm, min(stoplist)]
                ORFgff = [ALEchrm, ALE[1], 'ALE_ORF', str(ALEstart - 1), str(stoppos[1] + 2), ALE[5], ALEstrand, ALE[7], ALE[8]] #THIS SHOULD INCLUDE STOPCODON
                
            elif ALEstrand == '-':
                stoppos = [ALEchrm, max(stoplist)]
                ORFgff = [ALEchrm, ALE[1], 'ALE_ORF', str(stoppos[1]-3), str(ALEstop), ALE[5], ALEstrand, ALE[7], ALE[8]] #THIS ALSO INCLUDES STOPCODON

            ORFgffAll.append(ORFgff)

        elif len(stoplist) == 0: # if no stop codon found
            print "WARNING: no stop codon found for " + ALE[8]

        '''#This will print chromosome, length of ALE in coding region, length of ALE in UTR
        if len(stoplist) > 0:
            if ALEstrand == '+':
                print ALEchrm + '\t' + str(stoppos[1] - ALEstart) + '\t' + str(ALEstop - stoppos[1])
            elif ALEstrand == '-':
                print ALEchrm + '\t' + str(ALEstop - stoppos[1]) + '\t' + str(stoppos[1] - ALEstart)
        '''
        stoplist = [] #when done with ALE, reset lists
        stoppos = []
        ORFgff = []

    return ORFgffAll

def gfftosequence(gff, sequence_file):

    seq_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, 'fasta'))
    
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
            sequence = seq_dict[chrm].seq[start:stop].upper() #can slice differently here...this removes stop codons
        elif strand == '-':
            sequence = seq_dict[chrm].seq[start:stop].reverse_complement().upper()

        seqs[ID] = sequence

        GCs.append(float(GC(sequence)))

        counter +=1
        if counter <= 50 and counter % 10 == 0:
            print 'Retrieving sequence %i of %i' % (counter,len(gffentries))
        elif counter > 50 and counter % 50 == 0:
            print 'Retrieving sequence %i of %i' % (counter,len(gffentries))

    print 'The average GC content of these sequences is %f' % np.mean(GCs)

    return seqs

def translateALECDS(fasta):
    seqs = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        ntseq = record.seq
        ID = record.id
        if len(ntseq) % 3 == 0:
            aaseq = ntseq.translate()
        elif len(ntseq) % 3 == 1:
            ntseq = ntseq[1:]
            aaseq = ntseq.translate()
        elif len(ntseq) % 3 == 2:
            ntseq = ntseq[2:]
            aaseq = ntseq.translate()

        if aaseq.count('*') == 1: #if only one stop codon in translated seq:
            seqs[ID] = aaseq

    return seqs
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ALEgff', type = str, help = 'Gff of ALE events. Usually from MISO annotations.')
    parser.add_argument('--stopcodongff', type = str, help = 'Gff containing stop codons explicitly annotated. mm9.genes.chr.jess.gff works well.')
    parser.add_argument('--gffout', type = str, help = 'Outfile of gff of coding portion of ALE.')
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format.')
    parser.add_argument('--fastaout', type = str, help = 'Outfile of sequences of coding portion of ALE')
    parser.add_argument('--translatedout', help = 'Optional. If supplied, will translate CDS sequences and output in translatedout.')
    args = parser.parse_args()
       
    gfffh = open(args.gffout, 'w')
    print 'Determining coding portions of ALEs...'

    ALEORFCoords = getALEORFCoords(args.ALEgff, args.stopcodongff)

    for entry in ALEORFCoords:
        gfffh.write(('\t').join(entry) + '\n') #Output gff of coding portion of ALEs
        
    gfffh.close()

    print 'Retreiving sequences of coding portions of ALEs...'
    ORFseqs = gfftosequence(args.gffout, args.genomefasta)
    
    seqfh = open(args.fastaout, 'w')
    for ID in ORFseqs:
        seqfh.write('>' + ID + '\n' + str(ORFseqs[ID]) + '\n') #print in fasta format
    seqfh.close()

    if args.translatedout:
        aafh = open(args.translatedout, 'w')
        aaseqs = translateALECDS(args.fastaout)
        for ID in aaseqs:
            aafh.write('>' + ID + '\n' + str(aaseqs[ID]) + '\n') #print out aa seqs in fasta format
        aafh.close()
    
