#Takes fasta as input.  Calculates folding energies of 40 bp windows along sequence, with steps of 20 bp.  Reports median folding energy back as the energy for that sequence. Requires vienna RNA package (module load viennarna).

#Returns list of [ID, median_folding_energy].  As standalone, prints this list with tab separations.

#Necessary modules: biopython viennarna

#Usage: python UTRfold.py --help

import sys
from Bio import SeqIO
import fold_a_seq
import numpy as np
import argparse

def UTRfold(fasta):
    seqs = []
    medianenergies = []
    allenergies = []
    counter = 0

    for seq_record in SeqIO.parse(fasta, 'fasta'):
        record = []
        record.append(seq_record.id)
        record.append(str(seq_record.seq.transcribe()))
        seqs.append(record)

    for seq in seqs:
        currentlocation = 0
        ID = seq[0]
        sequence = seq[1]
        energies = []
        if len(sequence) >= 40:
            while currentlocation + 40 <= len(sequence):
                currentseq = sequence[currentlocation:currentlocation + 40]
                foldingenergy = float(fold_a_seq.fold_a_seq(currentseq))
                energies.append(foldingenergy)
                currentlocation += 20

            medianenergies.append([ID, np.median(energies)])
            allenergies.append([ID, energies])

        if counter % 100 == 0:
            print 'Folding sequence %i of %i' % (counter, len(seqs))
        counter +=1

    print 'Folded %i sequences.  %i were too short to fold.' % (len(medianenergies), len(seqs) - len(medianenergies))


    return allenergies


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to fold.')
    #parser.add_argument('--UTRtype', type = str, help = 'Type of UTRs in file (Localized distal, etc.)')
    parser.add_argument('--output', type = str, help = 'Output file.')
    args = parser.parse_args()
    medianenergies = UTRfold(args.fasta)
    outfh = open(args.output, 'w')
    for entry in medianenergies:
        #outfh.write(str(entry[0]) + '\t' + str(entry[1]) + '\t' + args.UTRtype + '\n')
        for energy in entry[1]:
            outfh.write(str(entry) + '\t' + str(energy) + '\n')
        
    outfh.close()
