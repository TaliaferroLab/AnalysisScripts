#Designed to take homologous regions of a genome (in this case 3' UTRs) and align those regions,
#then use RNAalifold to calculate folding energies for those alignments.  Takes as input bed files
#from mm9, rn5, hg19, canFam2, and bosTau7 as well as fasta files for their genomes.
#These beds were made by UCSC liftOver of coordinates of mm9 UTRs of interest.  Outputs the median energy
#of MFE of 100 nt windows of each alignment.

#Necessary modules: clustalw, biopython, viennarna

#Usage: python UTRfold_alignment.py --help

from Bio import SeqIO
from Bio import AlignIO
import sys
import subprocess
import numpy as np
import os
from Bio.Align.Applications import ClustalwCommandline
import replace_in_file
import argparse



def RNAalifold_a_seq(alignment):
    command = 'RNAalifold ' + alignment
    job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    output = job.communicate()
    return float(output[0].split('=')[0].split('(')[-1]) #This is the MFE from RNAalifold output

def foldAlignment(clustalaln):
    alignment = AlignIO.read(clustalaln, 'clustal')
    for record in alignment:
        if record.id[:3] == 'mm9':
            ID = record.id[4:]
    energies = []
    alnlength = alignment.get_alignment_length()
    currentlocation = 0
    print 'Folding {0}...'.format(ID)
    if alnlength >= 100:
        while currentlocation + 100 <= alnlength:
            #Take all rows (species) and only the columns (positions) of interest
            alignmentslice = alignment[:, currentlocation:currentlocation + 100]
            #Write the slice to a temporary file
            AlignIO.write(alignmentslice, 'tempslice.aln', 'clustal')
            
            #Align the slice, retrieving the MFE for the slice
            energy = RNAalifold_a_seq('tempslice.aln')
            energies.append(energy)
            
            #Move forward 100 nt
            currentlocation += 25
        os.remove('tempslice.aln')  
        return ID, np.median(energies)
    else:
        print 'Could not fold {0} because it was too short.'.format(ID)
        return ID, 'NA'
    
        
#May not actually need this function?
def getHomologousBeds(mm9bed, rn5bed, hg19bed, canFam2bed, bosTau7bed):
    mm9beds = []
    rn5beds = []
    hg19beds = []
    canFam2beds = []
    bosTau7beds = []
    nonmm9beds = [] #list of all non-mm9 beds (list of lists)
    bedhomologydict = {} #{ID: [{mm9 : [bed]}, {rn5 : [bed]}, {hg19 : [bed]}]}
    
    mm9bedfh = open(mm9bed, 'r')
    hg19bedfh = open(hg19bed, 'r')
    rn5bedfh = open(rn5bed, 'r')
    canFam2bedfh = open(canFam2bed, 'r')
    bosTau7bedfh = open(bosTau7bed, 'r')
    for line in mm9bedfh:
        line = line.strip().split('\t')
        mm9beds.append(line)
    for line in hg19bedfh:
        line = line.strip().split('\t')
        hg19beds.append(line)
    for line in rn5bedfh:
        line = line.strip().split('\t')
        rn5beds.append(line)
    for line in canFam2bedfh:
        line = line.strip().split('\t')
        canFam2beds.append(line)
    for line in bosTau7bedfh:
        line = line.strip().split('\t')
        bosTau7beds.append(line)
        
    mm9bedfh.close()
    hg19bedfh.close()
    rn5bedfh.close()
    canFam2bedfh.close()
    bosTau7bedfh.close()

    nonmm9beds = [rn5beds, hg19beds, canFam2beds, bosTau7beds]

    #For every bed in mm9, get homologous beds based on ID's (bed[3])
    for entry in mm9beds:
        mm9ID = entry[3][3:]
        mouseentry = {}
        mouseentry['mm9'] = entry
        homologousbeds = [] #[{mm9 : [bed]}, {rn5 : [bed]}, {hg19 : [bed]}]
        homologousbeds.append(mouseentry)
        for idx, species in enumerate(nonmm9beds):
            if idx == 0:
                genome = 'rn5'
            elif idx == 1:
                genome = 'hg19'
            elif idx == 2:
                genome = 'canFam2'
            elif idx == 3:
                genome = 'bosTau7'

            for bed in species:
                ID = bed[3][3:]
                if ID == mm9ID:
                    #This relies on there being one to one matches. There cannot be 2 matches in a species for 1 mouse entry.
                    speciesentry = {} 
                    speciesentry[genome] = bed
                    homologousbeds.append(speciesentry)

        bedhomologydict[mm9ID] = homologousbeds
    
    return bedhomologydict


#From dictionary of homologous UTR beds, get sequence for each bed.  Return in a similarly-structured dictionary.
def getSequences(bedhomologydict, mm9fasta, rn5fasta, hg19fasta, canFam2fasta, bosTau7fasta):
    sequencedict = {} #{UTR : [{'mm9_UTR' : sequence}, {'rn5_UTR' : sequence}, {'hg19_UTR' : sequence}]}
    
    sys.stderr.write('Indexing mm9 genome sequences...\n')
    mm9_seq_dict = SeqIO.to_dict(SeqIO.parse(mm9fasta, 'fasta'))
    sys.stderr.write('Indexed {0} mm9 sequences.\n'.format(len(mm9_seq_dict)))
    
    sys.stderr.write('Indexing rn5 genome sequences...\n')
    rn5_seq_dict = SeqIO.to_dict(SeqIO.parse(rn5fasta, 'fasta'))
    sys.stderr.write('Indexed {0} rn5 sequences.\n'.format(len(rn5_seq_dict)))

    sys.stderr.write('Indexing hg19 genome sequences...\n')
    hg19_seq_dict = SeqIO.to_dict(SeqIO.parse(hg19fasta, 'fasta'))
    sys.stderr.write('Indexed {0} hg19 sequences.\n'.format(len(hg19_seq_dict)))

    sys.stderr.write('Indexing canFam2 genome sequences...\n')
    canFam2_seq_dict = SeqIO.to_dict(SeqIO.parse(canFam2fasta, 'fasta'))
    sys.stderr.write('Indexed {0} canFam2 sequences.\n'.format(len(canFam2_seq_dict)))

    sys.stderr.write('Indexing bosTau7 genome sequences...\n')
    bosTau7_seq_dict = SeqIO.to_dict(SeqIO.parse(bosTau7fasta, 'fasta'))
    sys.stderr.write('Indexed {0} bosTau7 sequences.\n'.format(len(bosTau7_seq_dict)))
    
    for UTR in bedhomologydict:
        if len(bedhomologydict[UTR]) == 1: #if only one entry for a UTR (which means mouse UTR only, no homologs)
            continue
        elif len(bedhomologydict[UTR]) > 1:
            homologoussequences = [] # [{mm9_ID : mm9sequence}, {rn5_ID : rn5sequence}]
            for species in bedhomologydict[UTR]:
                if species.keys() == ['mm9']:
                    mouseentry = {}
                    chrm = species['mm9'][0]
                    start = int(species['mm9'][1])
                    stop = int(species['mm9'][2])
                    strand = species['mm9'][5]
                    mouseid = 'mm9_' + str(UTR)
                    if strand == '+':
                        UTRsequence = mm9_seq_dict[chrm].seq[start-1:stop].upper().transcribe()
                        mouseentry[mouseid] = str(UTRsequence)
                    elif strand == '-':
                        UTRsequence = mm9_seq_dict[chrm].seq[start-1:stop].upper().reverse_complement().transcribe()
                        mouseentry[mouseid] = str(UTRsequence)
                    homologoussequences.append(mouseentry)
                elif species.keys() == ['rn5']:
                    ratentry = {}
                    chrm = species['rn5'][0]
                    start = int(species['rn5'][1])
                    stop = int(species['rn5'][2])
                    strand = species['rn5'][5]
                    ratid = 'rn5_' + str(UTR)
                    if strand == '+':
                        UTRsequence = rn5_seq_dict[chrm].seq[start-1:stop].upper().transcribe()
                        ratentry[ratid] = str(UTRsequence)
                    elif strand == '-':
                        UTRsequence = rn5_seq_dict[chrm].seq[start-1:stop].upper().reverse_complement().transcribe()
                        ratentry[ratid] = str(UTRsequence)
                    homologoussequences.append(ratentry)
                elif species.keys() == ['hg19']:
                    humanentry = {}
                    chrm = species['hg19'][0]
                    start = int(species['hg19'][1])
                    stop = int(species['hg19'][2])
                    strand = species['hg19'][5]
                    humanid = 'hg19_' + str(UTR)
                    if strand == '+':
                        UTRsequence = hg19_seq_dict[chrm].seq[start-1:stop].upper().transcribe()
                        humanentry[humanid] = str(UTRsequence)
                    elif strand == '-':
                        UTRsequence = hg19_seq_dict[chrm].seq[start-1:stop].upper().reverse_complement().transcribe()
                        humanentry[humanid] = str(UTRsequence)
                    homologoussequences.append(humanentry)
                elif species.keys() == ['canFam2']:
                    dogentry = {}
                    chrm = species['canFam2'][0]
                    start = int(species['canFam2'][1])
                    stop = int(species['canFam2'][2])
                    strand = species['canFam2'][5]
                    dogid = 'canFam2_' + str(UTR)
                    if strand == '+':
                        UTRsequence = canFam2_seq_dict[chrm].seq[start-1:stop].upper().transcribe()
                        dogentry[dogid] = str(UTRsequence)
                    elif strand == '-':
                        UTRsequence = canFam2_seq_dict[chrm].seq[start-1:stop].upper().reverse_complement().transcribe()
                        dogentry[dogid] = str(UTRsequence)
                    homologoussequences.append(dogentry)
                elif species.keys() == ['bosTau7']:
                    cowentry = {}
                    chrm = species['bosTau7'][0]
                    start = int(species['bosTau7'][1])
                    stop = int(species['bosTau7'][2])
                    strand = species['bosTau7'][5]
                    cowid = 'bosTau7_' + str(UTR)
                    if strand == '+':
                        UTRsequence = bosTau7_seq_dict[chrm].seq[start-1:stop].upper().transcribe()
                        cowentry[cowid] = str(UTRsequence)
                    elif strand == '-':
                        UTRsequence = bosTau7_seq_dict[chrm].seq[start-1:stop].upper().reverse_complement().transcribe()
                        cowentry[cowid] = str(UTRsequence)
                    homologoussequences.append(cowentry)

            sequencedict[UTR] = homologoussequences #add homologous sequences to dictionary

    return sequencedict

def alignSeqs(sequencedict):
    energydict = {} #{UTRID : median_folding_energy}
    clustalfh = open('clustal_alignments.aln', 'w')
    UTRfastasfh = open('UTRfastas.fa', 'w')
    clustalfh.close()
    UTRfastasfh.close()
    counter = 0

    for UTR in sequencedict:
        counter +=1
        if counter % 10 == 0:
            print 'Folding sequence {0} of {1}...'.format(counter, len(sequencedict))
        UTRID = str(UTR)
        #Write fasta file from dictionary entry
        fastafh = open('temp.fasta', 'w')
        fastastring = ''
        for species in sequencedict[UTR]:
            fastastring += '>' + str(species.keys()[0]) + '\n' + str(species.values()[0]) + '\n'

        fastafh.write(fastastring)
        fastafh.close()
        tempfastafh = open('temp.fasta', 'r')
        tempfastalines = []
        for line in tempfastafh:
            tempfastalines.append(line)
        tempfastafh.close()

        #Align fasta using clustalw
        cline = ClustalwCommandline('clustalw2', infile = 'temp.fasta')
        cline() #alignment now in temp.aln
        clustallines = []
        tempclustalfh = open('temp.aln', 'r')
        for line in tempclustalfh:
            clustallines.append(line)
        tempclustalfh.close()

        #Fold temp.aln using foldAlignment. Returns median energy of 100 bp sections of alignment using RNAalifold
        ID, energy = foldAlignment('temp.aln')
        if energy:
            energydict[UTR] = energy

        #Append temp.aln to a big clustal file
        with open('clustal_alignments.aln', 'a') as clustalfile:
            for line in clustallines:
                clustalfile.write(line)

        #Append sequences to big fasta file of UTRs
        with open('UTRfastas.fa', 'a') as UTRfastafile:
            for line in tempfastalines:
                UTRfastafile.write(line)
            UTRfastafile.write('\n' + '\n' + '\n')

    #Cleanup
    os.remove('temp.aln')
    os.remove('temp.dnd')
    os.remove('temp.fasta')

    return energydict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mm9bed', type = str, help = 'Bed file of mm9 regions.')
    parser.add_argument('--rn5bed', type = str, help = 'Bed file of rn5 regions.')
    parser.add_argument('--hg19bed', type = str, help = 'Bed file of hg19 regions.')
    parser.add_argument('--canFam2bed', type = str, help = 'Bed file of canFam2 regions.')
    parser.add_argument('--bosTau7bed', type = str, help = 'Bed file of bosTau7 regions.')
    parser.add_argument('--mm9genome', type = str, help = 'mm9 genome in fasta format.')
    parser.add_argument('--rn5genome', type = str, help = 'rn5 genome in fasta format.')
    parser.add_argument('--hg19genome', type = str, help = 'hg19 genome in fasta format.')
    parser.add_argument('--canFam2genome', type = str, help = 'canFam2 genome in fasta format.')
    parser.add_argument('--bosTau7genome', type = str, help = 'bosTau7 genome in fasta format.')
    parser.add_argument('--output', type = str, help = 'Output file')
    args = parser.parse_args()

    bedhomologydict = getHomologousBeds(args.mm9bed, args.rn5bed, args.hg19bed, args.canFam2bed, args.bosTau7bed)
    sequencedict = getSequences(bedhomologydict, args.mm9genome, args.rn5genome, args.hg19genome, args.canFam2genome, args.bosTau7genome)
    energydict = alignSeqs(sequencedict)
    outfh = open(args.output, 'w')
    outfh.write('Sequence' + '\t' + 'Median_energy' + '\n')
    for ID in energydict:
        outfh.write(str(ID) + '\t' + str(energydict[ID]) + '\n')
    outfh.close()
    

