#Take a fasta of sequences you are interested in and count how many times kmers you are interested in occur.  
#Then correlate that with the deltapsi from the event that those sequences came from.

#python kmerdoseresponse.py --help

from Bio import SeqIO
import argparse

def getdeltaPSIs(deltaPSItable):
    deltapsis = {}
    tablefh = open(deltaPSItable, 'r')
    for event in tablefh:
        event = event.strip().split('\t')
        eventname = event[0]
        if eventname != 'Event':
            deltadeltapsi = float(event[8]) #MAY NEED TO CHANGE TO FIT COLUMN OF DELTADELTAPSIs
            deltapsis[eventname] = deltadeltapsi
    tablefh.close()

    return deltapsis

def getkmerstrengths(kmerstrengthsfile):
    kmerstrengths = {} #{kmer : strength}
    infh = open(kmerstrengthsfile, 'r')
    for line in infh:
        line = line.strip().split('\t')
        if line[0] != 'Kmer': #skip header
            kmerstrengths[line[0].upper()] = int(line[4])

    infh.close()

    return kmerstrengths

def countKmers(fasta, kmerstrengths, deltapsis, k):
    kmercounts = {} #{seqid : [sumofkmerscore, seqlength]}
    outdict = {} #{seqid : [sumofkmerscore, seqlength, deltapsi]}
    currentkmer = ''
    currentkmercount = 0
    k = int(k)
    counter = 0

    for seq_record in SeqIO.parse(fasta, 'fasta'):
        recordid = seq_record.id.split(';')[0][:-2] #MAY NEED TO SLICE BASED ON RELATIONSHIP OF ID IN FASTA TO ID IN DELTAPSI TABLE
        seq = str(seq_record.seq) #Used to .transcribe() here but RBNS kmers are DNA

        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            #Count homopolymeric stretches as only one instance of kmer, but allow nonoverlapping instances of
            #kmer in homopolymers
            if kmer in kmerstrengths:
                if kmercounts.has_key(recordid):
                    kmercounts[recordid][0] += kmerstrengths[kmer]
                else:
                    kmercounts[recordid] = [kmerstrengths[kmer]]
        if kmercounts.has_key(recordid):
            kmercounts[recordid].append(len(seq))
        else:
            kmercounts[recordid] = [0, len(seq)]
        counter +=1
        if counter % 100 == 0:
            print 'Analyzing sequence {0}...'.format(counter)

    for record in kmercounts:
        if record in deltapsis:
            outdict[record] = [kmercounts[record][0], kmercounts[record][1], deltapsis[record]]

    
    return outdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type = str, help = 'Fasta file of sequences to look through', required = True)
    parser.add_argument('--deltapsis', type = str, help = 'Delta psi table of events', required = True)
    parser.add_argument('--kmerstrengths', type = str, help = 'Table of kmer strengths from RBNS. 1st field = kmer. 5th field = kmerscore.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    args = parser.parse_args()

    deltapsis = getdeltaPSIs(args.deltapsis)
    kmerstrengths = getkmerstrengths(args.kmerstrengths)
    k = int(len(kmerstrengths.keys()[0]))
    kmercounts = countKmers(args.fasta, kmerstrengths, deltapsis, k)

    outfh = open(args.outfile, 'w')
    outfh.write('event' + '\t' + 'kmerscore' + '\t' + 'seqlength' + '\t'+ 'density' + '\t' + 'deltapsi' + '\t' + 'densitybin' +'\n')
    for event in kmercounts:
        kmercount = kmercounts[event][0]
        seqlength = float(kmercounts[event][1])
        density = kmercount / seqlength
        deltapsi = str(kmercounts[event][2])
        
        if density == 0:
            densitybin = 'none'
        elif density < 0.075 and density > 0:
            densitybin = 'verylow'
        elif density >= 0.075 and density < 0.15:
            densitybin = 'low'
        elif density >= 0.15 and density < 0.8:
            densitybin = 'medium'
        '''
        elif density >= 0.15 and density < 0.2:
            densitybin = 'high'
        elif density >= 0.2:
            densitybin = 'veryhigh'
        '''

    
        '''
        if density == 0:
            densitybin = 'none'
        elif density > 0 and density < 0.01:
            densitybin = 'low'
        elif density >= 0.01 and density < 0.02:
            densitybin = 'medium'
        elif density >= 0.02:
            densitybin = 'high'
        '''
        if seqlength >= 100:
            outfh.write(event + '\t' + str(kmercount) + '\t' + str(seqlength) + '\t' + str(density) + '\t' + deltapsi + '\t' + densitybin +'\n')

    outfh.close()
    

