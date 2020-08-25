#This is designed to take a fasta of oligo sequences designed by Jason
#for array-based oligo synthesis and retrieve their coordinates using BLAT.
#It requires that BLAT be in your path and that you have a directory of fasta
#sequences, one for each chromosome, named like chr1.fa.

import os
import subprocess
from Bio import SeqIO
import sys

#Usage: python fasta2gff.py <fasta> <Directory of Chromosome fasta sequences> <output>


def blatseq(fasta, chrm, seqdir):
    successfulalign = False
    for chrmseq in os.listdir(seqdir):
        if os.path.basename(chrmseq).split('.')[0] == chrm:
            database = os.path.join(seqdir, chrmseq)

    command = ['blat ' + database + ' ' + fasta + ' temp.psl']
    subprocess.call(command, shell=True)

    #Find the best match
    with open('temp.psl', 'r') as f:
        matchvalues = []
        for line in f:
            line = line.strip().split('\t')
            #Try and see if the first field is an integer. If it's not, try again on the next line.
            try:
                matchvalues.append(int(line[0]))
            except ValueError:
                continue

    if len(matchvalues) == 0: #if there is no alignment
        os.remove('temp.psl')
        print 'No alignment for sequence {0}!!'.format(os.path.basename(fasta))
        return None, None, None, successfulalign
        
    elif len(matchvalues) > 0:
        #The next time we go through the file and check if matchvalue == bestmatchvalue,
        #the matchvalues are going to be strings because they will have been split on tabs again,
        #so turn bestmatchvalue into a string to make matching easier
        bestmatchvalue = str(max(matchvalues))

    #Get coords of best match
    if bestmatchvalue:
        with open('temp.psl', 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] == bestmatchvalue:
                    mismatches = int(line[1])
                    if mismatches != 0: #if there are mismatches
                        print 'Error, sequence {0} aligned with mismatches!!'.format(infoline[9])
                        successfulalign = 'mismatch'
                    else:
                        successfulalign = True

                    strand = line[8]
                    start = int(line[15])
                    end = int(line[16])

        os.remove('temp.psl')
        
        return start, end, strand, successfulalign

def fasta2gff(fasta, seqdir, output):
    outfh = open(output, 'w')
    outfh.close()
    consideredseqs = 0
    mismatchaligns = 0
    failedaligns = 0
    successfulaligns = 0
    
    for record in SeqIO.parse(fasta, 'fasta'):
        consideredseqs +=1
        record = record[:] #slice record to cut off adapters
        SeqIO.write(record, 'temp.fasta', 'fasta')
        if 'chr' in record.id:
            chrm = record.id.split(':')[0]
        elif 'chr' not in record.id:
            #chrm = 'chr' + record.id.split(':')[0]
            chrm = 'chr' + record.id.split('|')[4]
            age = record.id.split('|')[0]
            location = record.id.split('|')[-1]
            #if record.id.split('|')[8].endswith('_exon'): #only consider flanking exon oligos
        start, end, strand, successfulalign = blatseq('temp.fasta', chrm, seqdir)
        
        if successfulalign == True:
            successfulaligns +=1
        elif successfulalign == False:
            failedaligns +=1
        elif successfulalign == 'mismatch':
            mismatchaligns +=1
            
        ID = 'ID=' + record.id

        if successfulalign:
            with open(output, 'a') as f:
                f.write(('\t').join([chrm, age, location, str(start), str(end), '.', strand, '.', ID]) + '\n')

        os.remove('temp.fasta')

    print 'Tried to align {0} sequences. {1} aligned successfully. {2} had at least one mismatch in the best alignment. {3} had no alignment at all.'.format(consideredseqs, successfulaligns, mismatchaligns, failedaligns)


fasta2gff(sys.argv[1], sys.argv[2], sys.argv[3])
