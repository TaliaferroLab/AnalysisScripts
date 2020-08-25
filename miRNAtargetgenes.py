#Takes a Summary_counts file from TargetScan and either a miRNA name or a seed sequence.  Returns list of that are targeted by that miRNA or seed sequence.  Optionally, takes a list of genes of interest and will make special note of those genes that are targeted.

#Usage: python miRNAtargetgenes.py --help

import sys
import os
import argparse

def geneswithsites(targetscanfile, genelist, queryseedsite, querymiRNA):
    genes_with_sites = []
    interesting_genes_with_sites = []
    interesting_genes_in_targetscan = 0
    genes_with_seedsites = []
    interesting_genes_with_seedsites = []
    genes_with_miRNAtargets = []
    interesting_genes_with_miRNAtargets = []
    targetscanfh = open(targetscanfile, 'r')
    targetscanlines = []
    genes_with_queryseedsite = []
    for line in targetscanfh:
        line = line.strip().split('\t')
        if line[0] != 'Transcript ID': #skip header line
            targetscanlines.append(line)
    targetscanfh.close()

    if genelist:
        genelistfh = open(genelist, 'r')
        genes_of_interest = []
        for line in genelistfh:
            line = line.strip()
            genes_of_interest.append(line)
        genelistfh.close()

        for gene in genes_of_interest:
            for line in targetscanlines:
                if gene == line[1]:
                    interesting_genes_in_targetscan +=1
                    break
        sys.stderr.write('Found {0} of {1} supplied genes in TargetScan file.\n'.format(interesting_genes_in_targetscan, len(genes_of_interest)))

    for line in targetscanlines:
        transcript = line[0]
        gene = line[1]
        seed = line[2]
        targetsites = int(line[4]) + int(line[5]) + int(line[6]) + int(line[7]) + int(line[8]) + int(line[9]) + int(line[10]) + int(line[11])
        miRNA = line[12]
        if queryseedsite:
            if queryseedsite == seed and targetsites > 0:
                genes_with_seedsites.append(gene)
                if genelist:
                    if gene in genes_of_interest:
                        interesting_genes_with_seedsites.append(gene)
        if querymiRNA:
            if querymiRNA == miRNA and targetsites > 0:
                genes_with_miRNAtargets.append(gene)
                if genelist:
                    if gene in genes_of_interest:
                        interesting_genes_with_miRNAtargets.append(gene)

    if queryseedsite:
        if genelist:
            #Remove duplicates. They could be hits on multiple transcripts of the same gene.
            interesting_genes_with_seedsites = list(set(interesting_genes_with_seedsites))
            genes_with_seedsites = list(set(genes_with_seedsites))
            return genes_with_seedsites, interesting_genes_with_seedsites
        else:
            genes_with_seedsites = list(set(genes_with_seedsites))
            return genes_with_seedsites, None # just easier to always return two things
    if querymiRNA:
        if genelist:
            interesting_genes_with_miRNAtargets = list(set(interesting_genes_with_miRNAtargets))
            genes_with_miRNAtargets = list(set(genes_with_miRNAtargets))
            return genes_with_miRNAtargets, interesting_genes_with_miRNAtargets
        else:
            genes_with_miRNAtargets = list(set(genes_with_miRNAtargets))
            return genes_with_miRNAtargets, None
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--targetscan', type = str, help = 'TargetScan summary counts file')
    parser.add_argument('--genelist', type = str, help = 'Optional. List of genes of particular interest. E.g. Localized genes')
    parser.add_argument('--output', type = str, help = 'Output file.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--queryseedsite', type = str, help = 'Seed site to search for. Corresponds to 3rd field of TargetScan file.')
    group.add_argument('--querymiRNA', type = str, help = 'miRNA to search for.  Corresponds to 13th field of TargetScan file.')
    args = parser.parse_args()
    
    outfh = open(args.output, 'w')
    
    genes_with_hits, interesting_genes_with_hits = geneswithsites(args.targetscan, args.genelist, args.queryseedsite, args.querymiRNA)
    if args.genelist:
        outfh.write('***All genes with hits***\n\n')
        for gene in genes_with_hits:
            outfh.write(gene + '\n')
        outfh.write('\n***Interesting genes with hits***\n\n')
        for gene in interesting_genes_with_hits:
            outfh.write(gene + '\n')
    elif args.genelist == None:
        for gene in genes_with_hits:
            outfh.write(gene + '\n')
    outfh.close()
