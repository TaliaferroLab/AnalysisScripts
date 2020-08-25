#Reorder a MISO gff file.  Miso takes the inclusion isoform as the first listed mRNA.  This reorders the file by start coord.  All other fields (ID, etc.) are unchanged.  Designed to work on AFE and ALE gff3s to sort by "most distal" isoform.  This means most upstream isoform for AFEs and most downstream isoform for ALEs.

#Usage python ReorderGff3.py infile outfile

import sys
from operator import itemgetter

infh = open(sys.argv[1], 'r')
outfh = open(sys.argv[2], 'w')

events = []
currentgene = []
currentmrnas = []
currentexons = []

#Make a list of lists of lists (called events).  The innermost list is each line in the GFF.  The next lists gathers up all exons or mrnas belonging to the same event.  The next list up gathers up all lines belonging to the same event (gene).  The outer list (events) gathers up all events.

for line in infh:
    line = line.strip()
    line = line.split('\t')
    if len(line) == 9: #sometimes the gff has a header line...skip it
        
        if line[2] == 'gene':
            if len(currentgene) == 0:
                currentgene.append(line)
            elif len(currentgene) > 0:
                currentgene.append(currentmrnas)
                currentgene.append(currentexons)
                events.append(currentgene)
                currentgene = []
                currentmrnas = []
                currentexons = []
                currentgene.append(line)
        elif line[2] == 'mRNA':
            currentmrnas.append(line)
        elif line[2] == 'exon':
            currentexons.append(line)

currentgene.append(currentmrnas) 
currentgene.append(currentexons)
events.append(currentgene)

#for event in events:
#event[0] = gene, event[1] = the mrnas, event[2] = the exons

for event in events:
    strand = event[0][6]
    outfh.write(('\t').join(event[0]) + '\n')

    mrnas = event[1]
    if strand == '+': #+ for AFEs - for ALEs
        mrnas.sort(key=itemgetter(3)) #sort by start coord
        for mrna in mrnas:
            outfh.write(('\t').join(mrna) + '\n')
    elif strand == '-': #- for AFEs + for ALEs
        mrnas.sort(key=itemgetter(4), reverse=True) #reverse sort by stop coord
        for mrna in mrnas:
            outfh.write(('\t').join(mrna) + '\n')

    exons = event[2]
    if strand == '+': #+ for AFE - for ALEs
        exons.sort(key=itemgetter(3))
        for exon in exons:
            outfh.write(('\t').join(exon) + '\n')
    elif strand == '-': #- for AFEs + for ALEs
        exons.sort(key=itemgetter(4), reverse=True)
        for exon in exons:
            outfh.write(('\t').join(exon) + '\n')

infh.close()
outfh.close()
