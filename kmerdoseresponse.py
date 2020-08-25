#Take a fasta of sequences you are interested in and count how many times kmers you are interested in occur.  
#Then correlate that with the deltapsi from the event that those sequences came from.

#python kmerdoseresponse.py --help

from Bio import SeqIO
import argparse
import sys
import random
from math import log

def getdeltaPSIs(deltaPSItable):
    deltapsis = {}
    tablefh = open(deltaPSItable, 'r')
    for event in tablefh:
        event = event.strip().split('\t')
        eventname = event[0]
        if event[5] == 'NA' or event[6] == 'NA' or event[0] == 'Gene':
        	continue
        if 'Inf' in event[5] or 'Inf' in event[6]:
        	continue
        #if float(event[1]) < 5 or float(event[2]) < 5 or float(event[3]) < 5 or float(event[4]) < 5 or float(event[5]) < 1:
        	#continue
        if eventname != 'Event' and eventname != 'Gene':
            deltadeltapsi = (float(event[2]) - float(event[3])) - (float(event[4]) - float(event[5])) #MAY NEED TO CHANGE TO FIT COLUMN OF DELTADELTAPSIs
            #deltadeltapsi = log(float(event[6]), 2) - log(float(event[5]), 2)
            deltapsis[eventname] = deltadeltapsi
    tablefh.close()

    return deltapsis

def countKmers(fasta, kmersofinterest, deltapsis, k):
    kmercounts = {} #{seqid : [number of kmer occurences, seqlength]}
    outdict = {} #{seqid : [number of kmer occurences, seqlength, deltapsi]}
    currentkmer = ''
    currentkmercount = 0
    k = int(k)
    counter = 0

    for seq_record in SeqIO.parse(fasta, 'fasta'):
        recordid = seq_record.id.split(';')[0][:-2] #MAY NEED TO SLICE BASED ON RELATIONSHIP OF ID IN FASTA TO ID IN DELTAPSI TABLE
        #recordid = seq_record.id
        seq = str(seq_record.seq.transcribe())

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            #Count homopolymeric stretches as only one instance of kmer, but allow nonoverlapping instances of
            #kmer in homopolymers
            if kmer in kmersofinterest:
                if kmercounts.has_key(recordid):
                    kmercounts[recordid][0] +=1
                else:
                    kmercounts[recordid] = [1]
        if kmercounts.has_key(recordid):
            kmercounts[recordid].append(len(seq))
        else:
            kmercounts[recordid] = [0, len(seq)]

    for record in kmercounts:
        if record in deltapsis:
            outdict[record] = [kmercounts[record][0], kmercounts[record][1], deltapsis[record]]

    return outdict

def makecontrolkmers(originalkmer):
    #Given a kmer, return all possible kmers with Ns between every base
    #For example: UGCU -> UNGNCNUN with all possible Ns
    #Afterwards, throw out any kmer with a different number of CpGs than starting kmer
    ctrlkmers = []
    currentctrlkmers = []
    alphabet = ['C','G','A','U']
    originalCpGcount = originalkmer.count('CG')

    for i in range(len(originalkmer)):
        if len(ctrlkmers) == 0:
            ctrlkmers.append(originalkmer[i])
            for kmer in ctrlkmers:
                for nt in alphabet:
                    currentctrlkmers.append(kmer + nt)

            ctrlkmers = currentctrlkmers
            currentctrlkmers = []

        elif len(ctrlkmers) > 0:
            for kmer in ctrlkmers:
                for nt in alphabet:
                    currentctrlkmers.append(kmer + originalkmer[i] + nt)

            ctrlkmers = currentctrlkmers
            currentctrlkmers = []

    ctrlkmers = list(set(ctrlkmers))
    for kmer in ctrlkmers:
        if kmer.count('CG') != originalCpGcount:
            ctrlkmers.remove(kmer)
        try:
            if originalkmer in kmer:
                ctrlkmers.remove(kmer)
        except:
            print originalkmer, kmer

    print 'After CpG matching, end up with {0} control kmers for {1}.'.format(len(ctrlkmers), originalkmer)
    return ctrlkmers

def makerandomkmers(kmer, iterations):
    iterations = int(iterations)
    sys.stderr.write('Creating {0} random kmers with matched GC content and CpG counts for kmer {1}...\n'.format(iterations, kmer))
    controlsequences = []
    kmer = kmer.upper()
    CpGcount = kmer.count('CG')
    GCcontent = kmer.count('C') + kmer.count('G')
    kmersize = len(kmer)
    while len(controlsequences) < iterations:
        randomseq = ''.join(random.choice('ACGU') for x in range(kmersize))
        randomCpGcount = randomseq.count('CG')
        randomGCcontent = randomseq.count('C') + randomseq.count('G')
        if randomCpGcount == CpGcount and randomGCcontent == GCcontent and randomseq != kmer:
            controlsequences.append(randomseq)
            
    sys.stderr.write('Done!\n')
    return controlsequences
            
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to look through', required = True)
    parser.add_argument('--deltapsis', type = str, help = 'Delta psi table of events', required = True)
    parser.add_argument('--kmers', type = str, help = 'Comma separated list of kmers you care about. Must all be same length.', required = True)
    parser.add_argument('--usectrlkmers', required = False, action = 'store_true', help = 'If given, instead of searching for kmers given by --kmers, use them to create control kmers containing Ns after every kmer nucleotide and search for those.')
    parser.add_argument('--userandomkmers', type = int, required = False, help = 'Mutually exclusive with usectrlkmers. Number of random C + G and CpG matched kmers to use as controls.')
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    parser.add_argument('--bin', type = str)
    args = parser.parse_args()

    if ',' in args.kmers:
        kmers = args.kmers.split(',')
    elif ',' not in args.kmers:
        kmers = [args.kmers]
    kmers = [kmer.upper() for kmer in kmers]
    
    k = len(kmers[0])

    
    #If we are wanting to use control kmers for the given kmers instead...
    if args.usectrlkmers == True:
        allcontrolkmers = []
        for kmer in kmers:
            controlkmers = makecontrolkmers(kmer)
            for controlkmer in controlkmers:
                allcontrolkmers.append(controlkmer)

        controlkmers = list(set(allcontrolkmers))

        print 'After collapsing, end up with {0} unique control kmers.'.format(len(controlkmers))
        kmers = controlkmers
        k = len(kmers[0])

    if args.userandomkmers:
        allcontrolkmers = []
        for kmer in kmers:
            controlkmers = makerandomkmers(kmer, args.userandomkmers)
            for controlkmer in controlkmers:
                allcontrolkmers.append(controlkmer)

        controlkmers = list(set(allcontrolkmers))
        controlkmers = controlkmers[:args.userandomkmers] #take first n of list of controlkmres
        kmers = controlkmers
        k = len(kmers[0])
        print kmers


    deltapsis = getdeltaPSIs(args.deltapsis)
    kmercounts = countKmers(args.fasta, kmers, deltapsis, k)

    outfh = open(args.outfile, 'w')
    outfh.write('event' + '\t' + 'kmercount' + '\t' + 'seqlength' + '\t'+ 'density' + '\t' + 'deltapsi' + '\t' + 'densitybin' + '\t' + 'Class' + '\n')
    for event in kmercounts:
        kmercount = kmercounts[event][0]
        seqlength = float(kmercounts[event][1])
        density = kmercount / seqlength
        deltapsi = str(kmercounts[event][2])
        
        
        if density == 0:
            densitybin = 1
        elif density < 0.0015 and density > 0:
            densitybin = 2
        elif density >= 0.0015 and density < 0.003:
            densitybin = 3
        elif density >= 0.003 and density < 0.004:
            densitybin = 4
        elif density >= 0.004 and density < 0.012:
            densitybin = 5
        elif density >= 0.012:
            densitybin = 6
        '''
        if density == 0:
            densitybin = '1'
        elif density > 0 and density < 0.003:
            densitybin = '2'
        elif density >= 0.003 and density < 0.006:
            densitybin = '3'
        elif density >= 0.009:
            densitybin = '4'
        '''
        
        
        ####################################CHECK THIS LINE!!!!!##################################
        outfh.write(event + '\t' + str(kmercount) + '\t' + str(seqlength) + '\t' + str(density / len(kmers)) + '\t' + deltapsi + '\t' + str(densitybin) + '\t' + args.bin +'\n')

    outfh.close()
    

