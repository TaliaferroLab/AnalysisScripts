#This will take two fasta files and output records for IDs that were in common between the two files. It's made for using MISO-based IDs of isoforms, e.g. '74351@uc007xnf.1@uc007xne.1.A'.  Here the '.A' distinguishes this isoform from its partner in the same gene, presumably '74351@uc007xnf.1@uc007xne.1.B'.

#Will also filter for records that pass the length filter in BOTH files. It also chops off the last <lengthfilter> nt of the sequence.  If you don't want a length filter, supply a length filter of 0.

#Usage : requireBothExons.py <fasta1> <fasta2> <lengthfilter> <fasta1output> <fasta2output>

import sys
from Bio import SeqIO

def requireBothExons(fasta1, fasta2, lengthfilter):

    lengthfilter = int(lengthfilter)
    fasta1genes = []
    fasta1dict = {}
    fasta2genes = []
    fasta2dict = {}
    fasta1out = {}
    fasta2out = {}

    fasta1genes = [[seq_record.id, len(seq_record)] for seq_record in SeqIO.parse(fasta1, 'fasta')] #every sequence name and length in fasta1
    fasta2genes = [[seq_record.id[:-2], len(seq_record)] for seq_record in SeqIO.parse(fasta2, 'fasta')] #every sequence name with isoform identifier (.X) removed and length in fasta2

    for gene in fasta1genes: 
        ID = gene[0]
        seqlength = gene[1]
        fasta1dict[ID] = [seqlength] #dictionary of name:length

    for gene in fasta2genes:
        ID = gene[0]
        seqlength = gene[1]
        fasta2dict[ID] = seqlength

    for entry in fasta1dict:
        choppedentry = entry[:-2]
        if choppedentry in fasta2dict:
            fasta1dict[entry].append(fasta2dict[choppedentry]) #dictionary of name: [length in fasta1, length in fasta2]
        else:
            fasta1dict[entry].append(int(0))

    for gene in fasta1dict:
        if fasta1dict[gene][0] > lengthfilter and fasta1dict[gene][1] > lengthfilter: #change legnth filters here
            for seq_record in SeqIO.parse(fasta1, 'fasta'):
                if gene == seq_record.id:
                    if lengthfilter > 0:
                        fasta1out[gene] = str(seq_record.seq[:-lengthfilter])
                    elif lengthfilter == 0:
                        fasta1out[gene] = str(seq_record.seq)

    print 'There are %i entries in Fasta 1 and %i entries in Fasta 2' % (len(fasta1genes), len(fasta2genes))
    print 'They have %i entries in common that pass length filter in both files' % len(fasta1out)
    print 'Building Fasta 1 using only common entries that pass length filter'
    
    fasta1genes = [] #reset lists and do the same thing for the second fasta file
    fasta2genes = []
    fasta1dict = {}
    fasta2dict = {}
    fasta1genes = [[seq_record.id[:-2], len(seq_record)] for seq_record in SeqIO.parse(fasta1, 'fasta')] #every sequence name with isoform identifier (.X) removed and length in fasta1
    fasta2genes = [[seq_record.id, len(seq_record)] for seq_record in SeqIO.parse(fasta2, 'fasta')] #every sequence name and length in fasta2

    for gene in fasta1genes:
        ID = gene[0]
        seqlength = gene[1]
        fasta1dict[ID] = seqlength #dictionary of name: length

    for gene in fasta2genes:
        ID = gene[0]
        seqlength = gene[1]
        fasta2dict[ID] = [seqlength]

    for entry in fasta2dict:
        choppedentry = entry[:-2]
        if choppedentry in fasta1dict:
            fasta2dict[entry].append(fasta1dict[choppedentry]) #dictionary of name: [length in fasta1, length in fasta2]
        else:
            fasta2dict[entry].append(int(0))

    print 'Building Fasta 2 using only common entries that pass length filter'

    for gene in fasta2dict:
        if fasta2dict[gene][0] > lengthfilter and fasta2dict[gene][1] > lengthfilter: #change length filters here
            for seq_record in SeqIO.parse(fasta2, 'fasta'):
                if gene == seq_record.id:
                    if lengthfilter > 0:
                        fasta2out[gene] = str(seq_record.seq[:-lengthfilter])
                    elif lengthfilter == 0:
                        fasta2out[gene] = str(seq_record.seq)

    return fasta1out, fasta2out #return tuple (of length 2) of dictionaries
            



fasta1outfh = open(sys.argv[4], 'w')
fasta2outfh = open(sys.argv[5], 'w')

fasta1out, fasta2out = requireBothExons(sys.argv[1], sys.argv[2], sys.argv[3])

for entry in fasta1out:
    fasta1outfh.write('>' + entry + '\n' + fasta1out[entry] + '\n')

for entry in fasta2out:
    fasta2outfh.write('>' + entry + '\n' + fasta2out[entry] + '\n')


fasta1outfh.close()
fasta2outfh.close()
