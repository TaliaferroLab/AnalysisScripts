import os
import argparse
from Bio import SeqIO
import numpy as np
import sys
import operator

def getkmerpositions(kmerlist, seqs):
    #All kmers must be same length
    kmerpos = {} #{sequencename : position of kmer start}
    k = int(len(kmerlist[0]))

    print 'Finding kmer positions in RNA pool...'
    for seq_record in SeqIO.parse(seqs, 'fasta'):
        seqname = seq_record.id
        seq = str(seq_record.seq.transcribe())

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer in kmerlist:
                if seqname not in kmerpos:
                    kmerpos[seqname] = [i + 1]
                elif seqname in kmerpos:
                    kmerpos[seqname].append(i + 1)

    return kmerpos

def mappingfilter(reactivities, filt):
    #Only consider RNAs with <filt> stops in the + channel
    print 'Filtering RNAs by stop counts...'
    filt = int(filt)
    filteredRNAs = []
    reactivitiesfh = open(reactivities, 'r')
    reactivitiesdict = {} # {sequence : total_treated_mods}
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence': #skip header
            sequence = line[0]
            treatedmods = int(line[4])
            if sequence in reactivitiesdict:
                reactivitiesdict[sequence] += treatedmods
            elif sequence not in reactivitiesdict:
                reactivitiesdict[sequence] = treatedmods
    for RNA in reactivitiesdict:
        if reactivitiesdict[RNA] >= filt:
            filteredRNAs.append(RNA)

    reactivitiesfh.close()

    print '{0} of {1} RNAs pass stop count filter.'.format(len(filteredRNAs), len(reactivitiesdict))

    return filteredRNAs

def getreactivities(kmerpos, filteredRNAs, reactivities, k):
    #Get reactivities at kmers of interest in RNAs that passed coverage filter
    reactivitydict = {} #{RNAid : {position : reactivity}}
    treatedmodsdict = {} #{RNAid : {position : treated_mods}}
    reactivityfh = open(reactivities, 'r')
    for line in reactivityfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence' and line[6] != '-': #skip header and position 0 of every sequence
            RNAid = line[0]
            if RNAid not in reactivitydict:
                reactivitydict[RNAid] = {}
            if RNAid not in treatedmodsdict:
                treatedmodsdict[RNAid] = {}
            position = int(line[2])
            treated_mods = int(line[4])
            #normalizedtheta = float(line[8])
            #reactivitydict[RNAid][position] = normalizedtheta
            bp_prob = float(line[10])
            reactivitydict[RNAid][position] = bp_prob
            treatedmodsdict[RNAid][position] = treated_mods
    reactivityfh.close()

    kmerreactivitydict = {} #{RNAid : {kmerstart : average reactivity}}
    for RNA in reactivitydict:
        if RNA in filteredRNAs and RNA in kmerpos: #if it passes both read count filter and it has at least one kmer instance in it
            kmerreactivitydict[RNA] = {}
            for kmerinstance in kmerpos[RNA]:
                ntreactivities = []
                kmertreatedmods = 0
                kmerstart = kmerinstance
                kmerend = kmerinstance + k
                for nt in range(kmerstart, kmerend):
                    treated_mods = treatedmodsdict[RNA][nt]
                    kmertreatedmods += treated_mods
                for nt in range(kmerstart, kmerend):
                    reactivity = reactivitydict[RNA][nt]
                    ntreactivities.append(reactivity)
                kmerreactivity = np.mean(ntreactivities)
                if kmertreatedmods >= 5: #Change this if you want to have a filter on at least <this many> stops in the + channel in the kmer
                    kmerreactivitydict[RNA][kmerstart] = kmerreactivity

    
    #Print all kmer instances sorted by mean reactivity (or bp prob according to line ~70)
    reacts = {}
    #print kmerreactivitydict
    for RNA in kmerreactivitydict:
        reacts[RNA] = 0
        for kmerstart in kmerreactivitydict[RNA]:
            if kmerreactivitydict[RNA][kmerstart] > reacts[RNA]:
                reacts[RNA] = kmerreactivitydict[RNA][kmerstart]

    #print sorted(reacts.items(), key = operator.itemgetter(1))
    return kmerreactivitydict

def getRvalues(rvalues):
    rfh = open(rvalues, 'r')
    rdict = {} # {seqname : r}
    for line in rfh:
        line = line.strip().split('\t')
        if line[0] != 'species':
            seqname = line[1]
            rvalue = float(line[2])
            rdict[seqname] = rvalue

    rfh.close()
    return rdict

def getreactivities_window(kmerpos, filteredRNAs, reactivities, rdict, k):
    #Get reactivities at kmers of interest in RNAs that passed coverage filter
    #kmerpos = {} #{sequencename : [positions of kmer start]}
    reactivitydict = {} #{RNAid : {position : reactivity}}
    treatedmodsdict = {} #{RNAid : {position : treated_mods}}
    reactivityfh = open(reactivities, 'r')
    for line in reactivityfh:
        line = line.strip().split('\t')
        if line[0] != 'sequence' and line[6] != '-': #skip header and position 0 of every sequence
            RNAid = line[0]
            if RNAid not in reactivitydict:
                reactivitydict[RNAid] = {}
            if RNAid not in treatedmodsdict:
                treatedmodsdict[RNAid] = {}
            position = int(line[1])
            treated_mods = int(line[4])
            #theta = float(line[7])
            #reactivitydict[RNAid][position] = theta
            bp_prob = float(line[10])
            reactivitydict[RNAid][position] = bp_prob
            treatedmodsdict[RNAid][position] = treated_mods
    reactivityfh.close()
    kmerreactivitydict = {} #{BNS : {age : {location : {binnumber : [bpprob1, bpprob2, etc...]}}}} #bin 1 is at kmerstart - windowsize
    kmerreactivitydict['bound'] = {}
    kmerreactivitydict['unbound'] = {}


    for RNA in reactivitydict:
        BNS = ''
        rvalue = rdict[RNA]
        if rvalue < 2: #CHANGE THESE IF NECESSARY!!!!!
            BNS = 'unbound'
        elif rvalue >= 2:
            BNS = 'bound'
        if RNA in filteredRNAs and RNA in kmerpos: #if it bassed both read count filter and it has at least one kmer instance in it
            age = RNA.split('|')[0]
            location = RNA.split('|')[-1]
            windowsize = 20 # + / - 10 bp on either side of kmer
            for kmerinstance in kmerpos[RNA]:
                #Only consider kmers that will have full windows, don't consider those whose window size will be off the RNA or will go into the RT primer binding site
                if int(kmerinstance) - windowsize > 0 and int(kmerinstance) + windowsize <= 112:
                    if kmerreactivitydict[BNS].has_key(age) == False:
                        kmerreactivitydict[BNS][age] = {}
                    if kmerreactivitydict[BNS][age].has_key(location) == False:
                        kmerreactivitydict[BNS][age][location] = {}
                    windowstart = int(kmerinstance) - windowsize
                    windowstop = int(kmerinstance) + k + windowsize
                    windowbin = 1 #start at first windowbin, which is windowsize upstream of the beginning of the kmer
                    for bp in range(windowstart, windowstop + 1):
                        if kmerreactivitydict[BNS][age][location].has_key(windowbin) == False:
                            kmerreactivitydict[BNS][age][location][windowbin] = [reactivitydict[RNA][bp]]
                        elif kmerreactivitydict[BNS][age][location].has_key(windowbin) == True:
                            kmerreactivitydict[BNS][age][location][windowbin].append(reactivitydict[RNA][bp])
                        windowbin += 1 

    return kmerreactivitydict




    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--kmerlist', type = str, help = 'Comma separated list of RNA kmers (Us not Ts) of interest. Must be same length.', required = True)
    parser.add_argument('--RNApool', type = str, help = 'RNA pool in fasta format.', required = True)
    parser.add_argument('--mappingfilter', type = int, help = 'RNA stop count filter. Only consider RNAs with at least this many stops in the + channel.',
                        required = True)
    parser.add_argument('--reactivities', type = str, help = 'reactivities.out from spats run.', required = True)
    parser.add_argument('--mode', choices = ['averageoverkmer', 'metaaroundkmer'], help = 'Do you want average reactivity over the kmer or a meta analysis around the kmer?')
    parser.add_argument('--rvalues', type = str, help = 'R values for oligos in reactivities.')
    parser.add_argument('--output', type = str, help = 'Output file.', required = True)

    args = parser.parse_args()

    if args.mode == 'averageoverkmer':
        kmerlist = args.kmerlist.split(',')
        k = int(len(kmerlist[0]))
        kmerpos = getkmerpositions(kmerlist, args.RNApool)
        filteredRNAs = mappingfilter(args.reactivities, args.mappingfilter)
        kmerreactivitydict = getreactivities(kmerpos, filteredRNAs, args.reactivities, k)
        rvalues = getRvalues(args.rvalues)

        outfh = open(args.output, 'w')
        outfh.write('RNA' + '\t' + 'class' + '\t' + 'location' + '\t' + 'kmerpos' + '\t' + 'kmerreactivity' + '\t' + 'Rvalue' + '\t' + 'BNS' + '\n')
        for RNA in kmerreactivitydict:
            age = RNA.split('|')[0]
            location = RNA.split('|')[7]
            rvalue = rvalues[RNA]
            BNS = None
            if rvalue >= 1.5:
                BNS = 'bound'
            elif rvalue < 1.5:
                BNS = 'unbound'
            for kmerpos in kmerreactivitydict[RNA]:
                reactivity = kmerreactivitydict[RNA][kmerpos]
                outfh.write(RNA + '\t' + age + '\t' + location + '\t' + str(kmerpos) + '\t' + str(reactivity) + '\t' + str(rvalue) + '\t' + BNS + '\n')
        outfh.close()

    elif args.mode == 'metaaroundkmer':
        def asint(s):
            try: 
                return int(s), ''
        
            except ValueError: 
                print s
                return sys.maxint, s


        kmerlist = args.kmerlist.split(',')
        k = int(len(kmerlist[0]))
        kmerpos = getkmerpositions(kmerlist, args.RNApool)
        filteredRNAs = mappingfilter(args.reactivities, args.mappingfilter)
        rvalues = getRvalues(args.rvalues)
        kmerreactivitydict = getreactivities_window(kmerpos, filteredRNAs, args.reactivities, rvalues, k) #{BNS : {age : {location : {binnumber : [bpprob1, bpprob2, etc...]}}})
        outfh = open(args.output, 'w')
        outfh.write(('\t').join(['BNS','Age','location','bin','meanbpprob', 'stdev', 'stderr']) + '\n')
        for BNS in kmerreactivitydict:
            for age in kmerreactivitydict[BNS]:
                for location in kmerreactivitydict[BNS][age]:
                    for binnumber in sorted(kmerreactivitydict[BNS][age][location], key = asint):
                        meanbpprob = round(np.mean(kmerreactivitydict[BNS][age][location][binnumber]), 4)
                        stdev = round(np.std(kmerreactivitydict[BNS][age][location][binnumber]), 4)
                        stderr = round(stdev / np.sqrt(len(kmerreactivitydict[BNS][age][location][binnumber])), 4)
                        outfh.write(('\t').join([BNS, age, location, str(binnumber), str(meanbpprob), str(stdev), str(stderr)]) + '\n')

        outfh.close()


