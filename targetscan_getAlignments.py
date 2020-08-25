#This script is designed to take beds of homologous 3' UTRs, get their sequences, and align
#them using clustalw.  The output is formatted for use with perl scripts associated with 
#Targetscan. Specifically, this output has three tab-separated fields: (1) UTR name,
#(2) taxonomy ID (species), and (3) aligned sequence.  Also outputs a fasta file of
#each 3'UTR for each species (those that are in the bed file) as UTRfastas.fa.  Homologous
#bed files were created by first converting the relevant mm9 gff of UTRs to a bed file using
#gfftobed.py.  The homologous regions to each bed entry was then retrieved using UCSC liftover.

#Currently setup to take mouse, rat, human, and dog sequences. (Fairly) easily modifiable to take more.

#Much of this script was based on getStructuralAlignments.py.  

#Usage: targetscan_getAlignments.py --help

import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import replace_in_file
import argparse

#May not actually need this function?
def getHomologousBeds(mm9bed, rn5bed, hg19bed, canFam2bed):
    mm9beds = []
    rn5beds = []
    hg19beds = []
    canFam2beds = []
    nonmm9beds = [] #list of all non-mm9 beds (list of lists)
    bedhomologydict = {} #{ID: [{mm9 : [bed]}, {rn5 : [bed]}, {hg19 : [bed]}]}
    
    mm9bedfh = open(mm9bed, 'r')
    hg19bedfh = open(hg19bed, 'r')
    rn5bedfh = open(rn5bed, 'r')
    canFam2bedfh = open(canFam2bed, 'r')
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
        
    mm9bedfh.close()
    hg19bedfh.close()
    rn5bedfh.close()
    canFam2bedfh.close()

    nonmm9beds = [rn5beds, hg19beds, canFam2beds]

    #For every bed in mm9, get homologous beds based on ID's (bed[3])
    for entry in mm9beds:
        mm9ID = entry[3][3:] #expecting 'ID=...'
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
def getSequences(bedhomologydict, mm9fasta, rn5fasta, hg19fasta, canFam2fasta):
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

            sequencedict[UTR] = homologoussequences #add homologous sequences to dictionary

    return sequencedict

def alignSeqs(sequencedict, outfile):
    taxonomydict = {'hg19' : '9606', 'mm9' : '10090', 'rn5' : '10116', 'canFam2' : '9615'}
    clustalfh = open(outfile, 'w')
    UTRfastasfh = open('UTRfastas.fa', 'w')
    clustalfh.close()
    UTRfastasfh.close()

    #Make temporary fasta for clustalw
    for UTR in sequencedict:
        UTRID = str(UTR)
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
        alignment = AlignIO.read('temp.aln', 'clustal')
        with open(outfile, 'a') as clustalfile:
            for record in alignment:
                species = record.id.split('_')[0]
                clustalfile.write(UTRID + '\t' + taxonomydict[species] + '\t' + str(record.seq) + '\n')

        #Add sequences to fastafile
        with open('UTRfastas.fa', 'a') as UTRfastafile:
            for line in tempfastalines:
                UTRfastafile.write(line)
            UTRfastafile.write('\n' + '\n' + '\n')

    os.remove('temp.fasta')
    os.remove('temp.aln')
    os.remove('temp.dnd')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mousebed', type = str, help = 'mouse sequences to consider in bed format.')
    parser.add_argument('--ratbed', type = str, help = 'rat homologs of mm9 sequences. possibly from UCSC liftover.')
    parser.add_argument('--humanbed', type = str, help = 'human homologs of mm9 sequences. possibly from UCSC liftover.')
    parser.add_argument('--dogbed', type = str, help = 'dog homologs of mm9 sequences. possibly from UCSC liftover.')
    parser.add_argument('--mousefasta', type = str, help = 'mouse genome in fasta format')
    parser.add_argument('--ratfasta', type = str, help = 'rat genome in fasta format')
    parser.add_argument('--humanfasta', type = str, help = 'human genome in fasta format')
    parser.add_argument('--dogfasta', type = str, help = 'dog genome in fasta format')
    parser.add_argument('--outfile', type = str, help = 'Output alignment file ready to use with TargetScan scripts.')
    args = parser.parse_args()

    
    bedhomologydict = getHomologousBeds(args.mousebed, args.ratbed, args.humanbed, args.dogbed)
    sequencedict = getSequences(bedhomologydict, args.mousefasta, args.ratfasta, args.humanfasta, args.dogfasta)
    alignSeqs(sequencedict, args.outfile)
