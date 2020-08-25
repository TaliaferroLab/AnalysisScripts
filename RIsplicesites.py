import sys
import argparse
from Bio import SeqIO

def getRI5ss(events, genomefasta):
    #Last 3 nt of exon, first 6 nt of intron
    fiveprimesscoords = {} #{eventname : [chrm, start, stop, strand]}
    fiveprimessseqs = {} #{eventname : sequence}

    eventsfh = open(events, 'r')
    for line in eventsfh:
        eventname = line.strip()
        line = line.strip().split('@')
        upstreamexon = line[0]
        upstreamexon = upstreamexon.split(':')
        chrm = upstreamexon[0]
        strand = upstreamexon[3]
        if strand == '+':
            upstreamexonstart = int(upstreamexon[1])
            upstreamexonstop = int(upstreamexon[2]) #last nt of exon
        elif strand == '-':
            upstreamexonstart = int(upstreamexon[2])
            upstreamexonstop = int(upstreamexon[1]) #last nt of exon

        fiveprimesscoords[eventname] = [chrm, upstreamexonstart, upstreamexonstop, strand]

    eventsfh.close()

    #Get sequences
    seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
    for event in fiveprimesscoords:
        chrm = fiveprimesscoords[event][0]
        upstreamexonend = fiveprimesscoords[event][2]
        strand = fiveprimesscoords[event][3]
        if strand == '+':
            ssstart = upstreamexonend - 3
            ssend = upstreamexonend + 6
            sequence = seq_dict[chrm].seq[ssstart:ssend].upper()
        elif strand == '-':
            ssstart = upstreamexonend - 7
            ssend = upstreamexonend + 2
            sequence = seq_dict[chrm].seq[ssstart:ssend].reverse_complement().upper()

        fiveprimessseqs[event] = str(sequence)

    return fiveprimessseqs 

def getRI3ss(events, genomefasta):
    #Last 20 nt of intron, first 3 nt of exon
    threeprimesscoords = {} #{eventname : [chrm, start, stop, strand]}
    threeprimessseqs = {} #{eventname : sequence}

    eventsfh = open(events, 'r')
    for line in eventsfh:
        eventname = line.strip()
        line = line.strip().split('@')
        downstreamexon = line[1]
        downstreamexon = downstreamexon.split(':')
        chrm = downstreamexon[0]
        strand = downstreamexon[3]
        if strand == '+':
            downstreamexonstart = int(downstreamexon[1])
            downstreamexonstop = int(downstreamexon[2]) #last nt of exon
        elif strand == '-':
            downstreamexonstart = int(downstreamexon[2])
            downstreamexonstop = int(downstreamexon[1]) #last nt of exon

        threeprimesscoords[eventname] = [chrm, downstreamexonstart, downstreamexonstop, strand]

    eventsfh.close()

    #Get sequences
    seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
    for event in threeprimesscoords:
        chrm = threeprimesscoords[event][0]
        downstreamexonstart = threeprimesscoords[event][1]
        strand = threeprimesscoords[event][3]
        if strand == '+':
            ssstart = downstreamexonstart - 21
            ssend = downstreamexonstart + 2
            sequence = seq_dict[chrm].seq[ssstart:ssend].upper()
        elif strand == '-':
            ssstart = downstreamexonstart - 3
            ssend = downstreamexonstart + 20
            sequence = seq_dict[chrm].seq[ssstart:ssend].reverse_complement().upper()

        threeprimessseqs[event] = str(sequence)

    return threeprimessseqs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type = str, help = 'File of RI events to consider. Required.', required = True)
    parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format. Required.', required = True)
    parser.add_argument('--splicesite', type = int, help = 'Which splice site? 5 or 3? Required.', choices = [5,3], required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    args = parser.parse_args()

    if args.splicesite == 5:
        seqs = getRI5ss(args.events, args.genomefasta)
    elif args.splicesite == 3:
        seqs = getRI3ss(args.events, args.genomefasta)

    outfh = open(args.outfile, 'w')
    for entry in seqs:
        outfh.write('>' + entry + '\n' + seqs[entry] + '\n')
    outfh.close()
