#Designed to take homologous regions of a genome (in this case 3' UTRs) and perform structural alignments.
#Takes as input bed files from mm9, rn5, hg19, canFam2, and bosTau7 as well as fasta files for their genomes. 
#These beds were made by UCSC liftOver of coordinates of mm9 UTRs of interest.  
#The mm9_bed UTRs had already had their stop codons and final 50 nt removed. 
#This script WILL NOT remove stop codons or last 50 nt.
#Outputs:
# (1) Fasta files for each UTR in each species (UTRfastas.fa)
# (2) Clustal alignments for each UTR (clustal_alignments.aln)
# (3) Stockholm alignments with SS_cons line for input to Infernal (stockholm_alignments.aln)

#Necessary modules: clustalw, biopython, viennarna

#Usage: python getStructuralAlignments.py <mm9bed> <rn5bed> <hg19bed> <canFam2bed> <bosTau7bed> <mm9genome.fa> <rn5genome.fa>
# <hg19genome.fa> <canFam2genome.fa> <bosTau7genome.fa>

import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import replace_in_file

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

#Take the fasta-like dictionary of UTR sequences and align them using clustalw.
#Then take that alignment and get consensus RNA structure using RNAalifold.
def alignSeqs(sequencedict):
    clustalfh = open('clustal_alignments.aln', 'w')
    '''
    stockholmfh = open('stockholm_alignments.aln', 'w')
    '''
    UTRfastasfh = open('UTRfastas.fa', 'w')
    clustalfh.close()
    '''
    stockholmfh.close()
    '''
    UTRfastasfh.close()

    if os.path.exists('./StockholmAlignments/') == False:
        os.mkdir('./StockholmAlignments')
    
    for UTR in sequencedict:
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
        
        #Convert clustal to stockholm
        AlignIO.convert('temp.aln', 'clustal', 'tempstockholm.aln', 'stockholm') 

        #Get secondary structure line from RNAalifold
        ss = subprocess.check_output(['RNAalifold', 'temp.aln']).replace(' ', '\n', 1).split('\n')[-3]
        ssline = '#=GC SS_cons ' + ss + '\n' + '//' + '\n'
        
        #Replace '//' in stockholm file with secondary structure line
        replace_in_file.replace('tempstockholm.aln', '//', ssline)

        #Add ID line to file.  This is necessary for Infernal.
        titleline = '# STOCKHOLM 1.0'
        IDline = '#=GF ID ' + UTRID
        replacement = titleline + '\n' + IDline + '\n'
        replace_in_file.replace('tempstockholm.aln', '# STOCKHOLM 1.0', replacement)

        #Rename stockholm file
        os.rename('tempstockholm.aln', './StockholmAlignments/' + UTRID + '.aln')

        #Now making many small stockholm files instead of one big one.
        '''
        tempstockholmfh = open('tempstockholm.aln', 'r')
        stockholmlines = []
        for line in tempstockholmfh:
            stockholmlines.append(line)
        tempstockholmfh.close()
        '''
        
        #Append current temp aln files to their respective alignment files
        with open('clustal_alignments.aln', 'a') as clustalfile:
            for line in clustallines:
                clustalfile.write(line)

        #Now making many small stockholm files instead of one big one.
        '''
        with open('stockholm_alignments.aln', 'a') as stockholmfile:
            for line in stockholmlines:
                stockholmfile.write(line)
        '''
        
        with open('UTRfastas.fa', 'a') as UTRfastafile:
            for line in tempfastalines:
                UTRfastafile.write(line)
            UTRfastafile.write('\n' + '\n' + '\n')

    #Cleanup
    os.remove('alirna.ps')
    os.remove('temp.aln')
    os.remove('temp.dnd')
    os.remove('temp.fasta')
        
    
if __name__ == '__main__':
    bedhomologydict = getHomologousBeds(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    sequencedict = getSequences(bedhomologydict, sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    alignSeqs(sequencedict)
    
