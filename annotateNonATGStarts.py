#Takes a gff.  Entries annotated as 'start_codon' are checked against fasta files of genome sequence to see if they are 'ATG'.  If so, they are unchanged.  If not, their category is changed to 'nonATG_start_codon'.  Gives essentially the same gff file back. Requires a single fasta of the entire genome.

#~~~~!!!!Reads the entire genome into memory!!! Best run on coyote!!!!~~~~~~

#Usage: python annotateNonATGstarts.py <in_gfffile> <genomefasta> <out_gfffile>

import os
import sys
from Bio import SeqIO

def annotateNonATGStarts(gff, genomefasta, outfile):

    seqdict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
    print '{0} sequences in fasta file'.format(len(seqdict))

    gfffh = open(gff, 'r')
    outfh = open(outfile, 'w')

    number_of_starts = 0
    number_of_ATG_starts = 0
    number_of_nonATG_starts = 0

    for line in gfffh:
        line = line.strip().split('\t')
        chrm = line[0]
        cat = line[2]
        start = int(line[3])
        stop = int(line[4])
        strand = line[6]
        sequence = ''
        
        if cat != 'start_codon' or 'random' in chrm: #skip chr*_random entries
            outfh.write(('\t').join(line) + '\n')
        elif cat == 'start_codon' and 'random' not in chrm:
            number_of_starts +=1
            if strand == '+':
                sequence = str(seqdict[chrm].seq[start-1:stop].upper())
            elif strand == '-':
                sequence = str(seqdict[chrm].seq[start-1:stop].upper().reverse_complement())

            if sequence == 'ATG':
                number_of_ATG_starts +=1
                outfh.write(('\t').join(line) + '\n')
            elif sequence != 'ATG':
                number_of_nonATG_starts +=1
                line = [entry.replace('start_codon', 'nonATG_start_codon') for entry in line]
                outfh.write(('\t').join(line) + '\n')

            

            if number_of_starts > 0 and number_of_starts <= 1000 and number_of_starts % 100 == 0:
                print 'Analyzed {0} start codons. So far {1} were ATG and {2} were non-ATG.'.format(number_of_starts, number_of_ATG_starts, number_of_nonATG_starts)
            elif number_of_starts > 1000 and number_of_starts % 1000 == 0:
                print 'Analyzed {0} start codons. So far {1} were ATG and {2} were non-ATG.'.format(number_of_starts, number_of_ATG_starts, number_of_nonATG_starts)
            
    

    print 'Found {0} annotated start codons.  Of these, {1} were ATG and {2} were non-ATG.'.format(number_of_starts, number_of_ATG_starts, number_of_nonATG_starts)

    gfffh.close()
    outfh.close()
    

if __name__ == '__main__':
    annotateNonATGStarts(sys.argv[1], sys.argv[2], sys.argv[3])
            

    
