#Takes a list of polyA sites from Derti et al, e.g. GSM747481_mouse_brain.sites.clustered.bed
#Also takes a gff genome annotation that has stop codons explicitly annotated as such (stop_codon)
#First finds polyA sites that are within 10 kb of a stop codon (not gonna deal with UTRs longer than that)
#Then finds the stop codon that is closest to that polyA site such that you get the shortest possible UTR.
#Then annotates the UTRs as such and filters for events that only have at least 2 isoforms (TandemUTRs).
#This output file should then be ready for indexing by MISO.

import sys
import os
import gffutils
import argparse
from operator import itemgetter

def getpolyAsites(bed):
    #Should be bed of sites, for exapmle from GSM747481
    polyAsites = [] #[chrm, coord, strand]
    bedfh = open(bed, 'r')
    for line in bedfh:
        line = line.strip().split('\t')
        chrm = line[0]
        strand = line[5]
        if strand == '+':
            coord = line[1]
        elif strand == '-':
            coord = line[2]
        polyAsites.append([chrm, coord, strand])

    bedfh.close()

    return polyAsites

def getstopcodons(gff):
    stopcodoncoords = [] #[chrm, coord, strand, ID]
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    print 'Indexing gff...\n'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    stopcodons = db.features_of_type('stop_codon')
    numberofstopcodons = sum(1 for _ in stopcodons)
    print 'Found {0} stop codons in gff.\n'.format(numberofstopcodons)
    #Remake generator
    stopcodons = db.features_of_type('stop_codon')

    counter = 0
    for stopcodon in stopcodons:
        counter +=1
        if counter % 1000 == 0:
            print 'Indexing stop coord {0} of {1}.\n'.format(counter, numberofstopcodons)
        chrm = stopcodon.chrom
        strand = stopcodon.strand
        if strand == '+':
            coord = int(stopcodon.stop) + 1
        elif strand == '-':
            coord = int(stopcodon.start) - 1
        for parent in db.parents(stopcodon, featuretype = 'mRNA'):
            #this recreates the ID and Name entries for the mRNA that is the parent of the stop codon
            ID = parent.id + ';Name=' + parent.attributes['Name'][0] 
        stopcodoncoords.append([str(chrm), coord, str(strand), str(ID)])

        
    return stopcodoncoords

def intersectsites(polyAsites, stopcodons):
    #Get polyA sites that lie within 10kb of a stop codon
    print 'OK, we have {0} polyA sites and {1} stop codons. Intersecting...\n'.format(len(polyAsites), len(stopcodons))
    potentialUTRs = [] #[chrm, strand, stopcoord, Acoord, ID]
    for polyAsite in polyAsites:
        Achrm = polyAsite[0]
        Acoord = int(polyAsite[1])
        Astrand = polyAsite[2]
        if Astrand == '+':
            for stopcodon in stopcodons:
                if stopcodon[2] == '+':
                    stopchrm = stopcodon[0]
                    stopcoord = int(stopcodon[1])
                    stopstrand = stopcodon[2]
                    stopID = stopcodon[3]
                    if Achrm == stopchrm and Acoord - stopcoord <= 10000 and Acoord - stopcoord > 0:
                        #polyAsite must be within 10 kb of stop coord (max UTR length is 10 kb)
                        potentialUTRs.append([Achrm, Astrand, stopcoord, Acoord, stopID])
        elif Astrand == '-':
            for stopcodon in stopcodons:
                if stopcodon[2] == '-':
                    stopchrm = stopcodon[0]
                    stopcoord = int(stopcodon[1])
                    stopstrand = stopcodon[2]
                    stopID = stopcodon[3]
                    if Achrm == stopchrm and stopcoord - Acoord <= 10000 and stopcoord - Acoord > 0:
                        potentialUTRs.append([Achrm, Astrand, stopcoord, Acoord, stopID])
                        
    print 'Found {0} potential UTRs.'.format(len(potentialUTRs))
    return potentialUTRs

def collapseUTRs(potentialUTRs):
    #If a particular polyAcoord has multiple stopcoords associated with it, take the stopcoord closest to the polyAcoord (shortest UTR)
    UTRs = {} #{chrm;strand;Acoord;stopid : [stopcoords]}
    collapsedUTRs = []
    for UTR in potentialUTRs:
        chrm = UTR[0]
        strand = UTR[1]
        stopcoord = UTR[2]
        Acoord = UTR[3]
        stopID = UTR[4]
        UTRID = (';').join([chrm, strand, str(Acoord), stopID])
        if UTRID not in UTRs:
            UTRs[UTRID] = [stopcoord]
        elif UTRID in UTRs:
            UTRs[UTRID].append(stopcoord)

    for UTR in UTRs:
        if len(UTRs[UTR]) > 1: #if more than one stop coord for that Acoord
            stopcoords = UTRs[UTR]
            if ';+;' in UTR:
                stopcoord = max(stopcoords)
            elif ';-;' in UTR:
                stopcoord = min(stopcoords)
        elif len(UTRs[UTR]) == 1:
            stopcoord = UTRs[UTR][0]
        chrm = UTR.split(';')[0]
        strand = UTR.split(';')[1]
        Acoord = int(UTR.split(';')[2])
        stopID = UTR.split(';')[3]
        collapsedUTRs.append([chrm, strand, stopcoord, Acoord, stopID])

    print 'After collapsing, we are left with {0} UTRs.'.format(len(collapsedUTRs))
    return collapsedUTRs

def formatforgff(collapsedUTRs):
    #To make 'gene' entries in gff
    #Need longest UTR for each ID
    plusgenecoords = {} #{stopcoord : [polyAsite1, polyAsite2]}
    minusgenecoords = {}
    gfflines = [] #list of lists, each entry is one line in the gff
    minusUTRs = 0
    for UTR in collapsedUTRs:
        chrm = UTR[0]
        strand = UTR[1]
        stopcoord = UTR[2]
        polyAcoord = UTR[3]
        stopID = UTR[4]
        UTRID = (';').join([chrm, strand, str(stopcoord), stopID])
        if strand == '+':
            if UTRID not in plusgenecoords:
                plusgenecoords[UTRID] = polyAcoord
            elif UTRID in plusgenecoords:
                if strand == '+' and polyAcoord > plusgenecoords[UTRID]:
                    plusgenecoords[UTRID] = polyAcoord
        elif strand == '-':
            if UTRID not in minusgenecoords:
                minusgenecoords[UTRID] = polyAcoord
            elif UTRID in minusgenecoords:
                if strand == '-' and polyAcoord < minusgenecoords[UTRID]:
                    minusgenecoords[UTRID] = polyAcoord

    #Take gene entries, and then any UTRs that have a stopcoord in common with either the genestop
    #(for those on the minus strand) or a genestart (for those on the plus strand)
    for gene in plusgenecoords:
        genechrm = gene.split(';')[0]
        genestrand = gene.split(';')[1]
        genestart = int(gene.split(';')[2])
        genestop = int(plusgenecoords[gene])
        geneID = gene.split(';')[3]
        if genestrand == '+':
            geneline = [genechrm, 'TandemUTR', 'gene', str(genestart), str(genestop), '.', genestrand, '.', 'ID='+ geneID]
        gfflines.append(geneline)
        isoformlines = []
        exonlines = []
        isoformcounter = 0
        for UTR in collapsedUTRs:
            if UTR[1] == '+' and genestrand == '+': #considering + strand only
                chrm = UTR[0]
                strand = UTR[1]
                stopcoord = UTR[2]
                polyAcoord = UTR[3]
                if chrm == genechrm and int(stopcoord) == int(genestart):
                    isoformcounter +=1
                    isoformline = [chrm, 'TandemUTR', 'mRNA', str(stopcoord), str(polyAcoord), '.', strand, '.', 'ID='+ geneID + '.' + str(isoformcounter)
                                   + ';Parent=' + geneID]
                    if any(str(polyAcoord) in entry for entry in isoformlines) == False: #must be a new polyAsite
                        isoformlines.append(isoformline)
                    exonline = [chrm, 'TandemUTR', 'exon', str(stopcoord), str(polyAcoord), '.', strand, '.', 'ID='+ geneID + '.' + str(isoformcounter)
                                   + '.0;Parent=' + geneID + '.' + str(isoformcounter)]
                    if any(str(polyAcoord) in entry for entry in exonlines) == False:
                        exonlines.append(exonline)

        #Sort these isoforms and exons so that the longest isoform (or exon) comes first in the list
        for isoformline in sorted(isoformlines, key=itemgetter(4), reverse=True):
            gfflines.append(isoformline)
        for exonline in sorted(exonlines, key=itemgetter(4), reverse=True):
            gfflines.append(exonline)
    
    for gene in minusgenecoords:
        genechrm = gene.split(';')[0]
        genestrand = gene.split(';')[1]
        genestart = int(minusgenecoords[gene])
        genestop = gene.split(';')[2]
        geneID = gene.split(';')[3]
        if genestrand == '-':
            geneline = [genechrm, 'TandemUTR', 'gene', str(genestart), str(genestop), '.', genestrand, '.', 'ID='+ geneID]
        gfflines.append(geneline)
        isoformcounter = 0
        isoformlines = []
        exonlines = []
        for UTR in collapsedUTRs:
            if UTR[1] == '-' and genestrand == '-': #now considering minus strand
                chrm = UTR[0]
                strand = UTR[1]
                stopcoord = UTR[2]
                polyAcoord = UTR[3]
                if chrm == genechrm and int(genestop) == int(stopcoord):
                    isoformcounter +=1
                    isoformline = [chrm, 'TandemUTR', 'mRNA', str(polyAcoord), str(stopcoord), '.', strand, '.', 'ID='+ geneID + '.' + str(isoformcounter)
                                   + ';Parent=' + geneID]
                    if any(str(polyAcoord) in entry for entry in isoformlines) == False: #must be a new polyAsite
                        isoformlines.append(isoformline)
                    exonline = [chrm, 'TandemUTR', 'exon', str(polyAcoord), str(stopcoord), '.', strand, '.', 'ID='+ geneID + '.' + str(isoformcounter)
                                   + '.0;Parent=' + geneID + '.' + str(isoformcounter)]
                    if any(str(polyAcoord) in entry for entry in exonlines) == False:
                        exonlines.append(exonline)
        for isoformline in sorted(isoformlines, key=itemgetter(3), reverse=False):
            gfflines.append(isoformline)
        for exonline in sorted(exonlines, key=itemgetter(3), reverse=False):
            gfflines.append(exonline)

    return gfflines

def filtergfflines(gfflines):
    #filter so that only genes with at least two isoforms pass
    filteredgfflines = []
    currentotherlines = []
    for line in gfflines:
        if line[2] == 'gene':
            if len(currentotherlines) > 2:
                filteredgfflines.append(currentgeneline)
                for entry in currentotherlines:
                    filteredgfflines.append(entry)
            currentotherlines = []
            currentgeneline = line
        elif line[2] != 'gene':
            currentotherlines.append(line)

    numberofevents = 0
    for line in filteredgfflines:
        if line[2] == 'gene':
            numberofevents +=1
    print 'After filtering for events with at least two isoforms, we are left with {0} TandemUTR events.'.format(numberofevents)

    return filteredgfflines

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--polyAsites', type = str, help = 'File of polyA sites.', required = True)
    parser.add_argument('--stopcodons', type = str, help = 'Gff file with annotated stop codons.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output gff file containing gff of only tandem UTRs.', required = True)
    args = parser.parse_args()

    polyAsites = getpolyAsites(args.polyAsites)
    stopcodons = getstopcodons(args.stopcodons)
    potentialUTRs = intersectsites(polyAsites, stopcodons)
    collapsedUTRs = collapseUTRs(potentialUTRs)
    gfflines = formatforgff(collapsedUTRs)
    filteredgfflines = filtergfflines(gfflines)
    outfh = open(args.outfile, 'w')
    for line in filteredgfflines:
        line = ('\t').join(line)
        outfh.write(line + '\n')
    outfh.close()

    
