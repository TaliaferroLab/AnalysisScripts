#Runs RNALfold on sequences in fasta file.  For each seq, reutrns folding energy for most stable
#local structure as reported by RNALfold. These energies are then normalized by their length since
#longer structures get more chances to get a "max" value.

#Alternatively, will report number of regions that are more stable than a given threshold energy.
#This is also normalized by length.

#Alternatively, will report minimum free energy for entire sequence.  This is also normalized by length.

#Requires modules viennarna, biopython

#Usage: python RNALfold.py --help

import argparse
import subprocess
from Bio import SeqIO
import sys
import re

def RNALfold_moststablelocal(seq, span):
    #Takes an RNA sequence as input. Returns energy of most stable local structure.
    successfulfold = False
    command = ['RNALfold ' + '-L ' + str(span)]
    output = []
    energies = []
    job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    job.stdin.write(str(seq))
    
    for line in job.communicate():
        output.append(re.findall('\d+.\d+', line)) #find all numbers, append them to line
    try:
        output = output[0] #first item contains energies, second item is empty list for some reason
        for energy in output:
            if '.' in energy: #Energies are decimals. Indexes are integers.
                energies.append(float(energy))
        energies.pop() #Last item is always MFE for entire RNA
        lowestlocalenergy = '-' + str(max(energies)) #Minus sign had been taken off energies, re-adding it here
        successfulfold = True
    #if only energy was MFE, len(energies) is 1, after popping it there is 0
    except ValueError: 
        lowestlocalenergy = 'N/A'
    #for some reason some seqs do not return any energies, giving len(energies) of 0, which can't be popped
    except IndexError:  
        lowestlocalenergy = 'N/A'
    
    return lowestlocalenergy, successfulfold

def RNALfold_MFE(seq, span):
    #Takes an RNA sequence as input. Returns MFE of RNA as calculated by RNALfold.
    successfulfold = False
    command = ['RNALfold ' + '-L ' + str(span)]
    output = job = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE,
                                    stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    job.stdin.write(str(seq))
    try:
        output = job.communicate()[0].split('\n')[-2] #get final MFE in kcal/mol
        output = float(output.replace('(','').replace(')','')) #remove parentheses
        successfulfold = True
    except ValueError:
        output = 'N/A'
    except IndexError:
        output = 'N/A'

    return output, successfulfold

def RNALfold_threshold(seq, span, threshold):
    #Takes an RNA sequence as input. Returns number of local structures more stable than threshold.
    successfulfold = False
    threshold = abs(threshold) #Don't know whether user put in pos or neg threshold
    thresholdpassers = 0
    command = ['RNALfold ' + '-L ' + str(span)]
    output = []
    energies = []
    job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    job.stdin.write(str(seq))

    for line in job.communicate():
        output.append(re.findall('\d+.\d+', line)) #find all numbers, append them to line
    try:
        output = output[0] #first item contains energies, second item is empty list for some reason
        for energy in output:
            if '.' in energy: #Energies are decimals. Indexes are integers.
                energies.append(float(energy))
        energies.pop() # Last item is always MFE for entire RNA
        for energy in energies:
            if energy >= threshold: #Count number that pass threshold
                thresholdpassers +=1
        successfulfold = True
    except ValueError:
        print 'Value Error!!!'
        thresholdpassers = 'N/A'
    except IndexError:
        print 'Index Error!!!'
        thresholdpassers = 'N/A'

    return thresholdpassers, successfulfold

def foldlongseqs(fasta, span, threshold, mode):
    energydict = {}
    successfulfolds = 0
    totalseqs = 0
    for record in SeqIO.parse(fasta, 'fasta'):
        sys.stderr.write('Folding {0}...\n'.format(record.id))
        totalseqs +=1
        recordID = record.id
        seq = record.seq
        seqlength = len(seq)
        if mode == 'localenergy':
            lowestlocalenergy, successfulfold = RNALfold_moststablelocal(seq, span)
            energydict[recordID] = [lowestlocalenergy, seqlength]
        elif mode == 'threshold':
            thresholdpassers, successfulfold = RNALfold_threshold(seq, span, threshold)
            energydict[recordID] = [thresholdpassers, seqlength]
        elif mode == 'MFE':
            MFE, successfulfold = RNALfold_MFE(seq, span)
            energydict[recordID] = [MFE, seqlength]
        if successfulfold == True:
            successfulfolds +=1
    
    return energydict, totalseqs, successfulfolds
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required = True)
    parser.add_argument('--inputfasta', type = str, help = 'Fasta file of seqs to fold.')
    parser.add_argument('--span', type = int, help = 'Max span length for basepaired nucleotides. Option -L for RNALfold. Default is 150')
    parser.add_argument('--UTRtype', type = str, help = 'UTR class.')
    parser.add_argument('--output', type = str, help = 'Output file.')
    group.add_argument('--localenergy', action = 'store_true', help = 'Use to calculate the most stable local structure in seq. Energy will be normalized by seq length.')
    group.add_argument('--threshold', type = int, help = 'Use to count number of regions that have structure more stable than <threshold>. Number of passing regions will be normalized by seq length.')
    group.add_argument('--MFE', action = 'store_true', help = 'Use to calculate minimum free energy for entire sequence.  Energy will be normalized by seq length.')
    
    args = parser.parse_args()
    if args.span == None:
        span = 150 #Set default
    else:
        span = args.span
    if args.localenergy:
        energydict, totalseqs, successfulfolds = foldlongseqs(args.inputfasta, span, None, 'localenergy')
        outfh = open(args.output, 'w')
        outfh.write('Sequence' + '\t' + 'Stablest_Local_Structure' + '\t' + 'Length' + '\t' + 'Normalized_energy' + '\t' + 'UTRtype' +'\n')
        for seq in energydict:
            if energydict[seq][0] == 'N/A':
                outfh.write(str(seq) + '\t' + 'N/A' + '\t' + str(int(energydict[seq][1])) + '\t' + 'N/A' + '\n')
            elif energydict[seq][0] != 'N/A':
                energy = float(energydict[seq][0])
                seqlength = float(energydict[seq][1])
                normalizedenergy = energy / seqlength
                outfh.write(str(seq) + '\t' + str(energy) + '\t' + str(int(seqlength)) + '\t' + str(normalizedenergy) + '\t' + args.UTRtype + '\n')

    elif args.threshold:
        energydict, totalseqs, successfulfolds = foldlongseqs(args.inputfasta, span, args.threshold, 'threshold')
        outfh = open(args.output, 'w')
        outfh.write('Sequence' + '\t' + 'Regions_that_pass_threshold' + '\t' + 'Length' + '\t' + 'Normalized_regions' + '\t' + 'UTRtype' + '\n')
        for seq in energydict:
            if energydict[seq][0] == 'N/A':
                outfh.write(str(seq) + '\t' + 'N/A' + '\t' + str(int(energydict[seq][1])) + '\t' + 'N/A' + '\n')
            elif energydict[seq][0] != 'N/A':
                regions = float(energydict[seq][0])
                seqlength = float(energydict[seq][1])
                normalizedregions = regions / seqlength
                outfh.write(str(seq) + '\t' + str(int(regions)) + '\t' + str(int(seqlength)) + '\t' + str(normalizedregions) +  '\t' + args.UTRtype + '\n')

    elif args.MFE:
        energydict, totalseqs, successfulfolds = foldlongseqs(args.inputfasta, span, None, 'MFE')
        outfh = open(args.output, 'w')
        outfh.write('Sequence' + '\t' + 'MFE' + '\t' + 'Length' + '\t' + 'Normalized_MFE' + '\t' + 'UTRtype' + '\n')
        for seq in energydict:
            if energydict[seq][0] == 'N/A':
                outfh.write(str(seq) + '\t' + 'N/A' + '\t' + str(int(energydict[seq][1])) + '\t' + 'N/A' + '\n')
            elif energydict[seq][0] != 'N/A':
                MFE = float(energydict[seq][0])
                seqlength = float(energydict[seq][1])
                normalizedMFE = MFE / seqlength
                outfh.write(str(seq) + '\t' + str(MFE) + '\t' + str(int(seqlength)) + '\t' + str(normalizedMFE) + '\t' + args.UTRtype + '\n')
            
    outfh.close()
    sys.stderr.write('Successfully folded {0} of {1} sequences.\n'.format(successfulfolds, totalseqs))
