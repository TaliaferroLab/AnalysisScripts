import os
import argparse
from math import log

def getFPKMratios(fpkm_tracking, comparisons, fpkmcutoff):
    FPKMs = {} # {gene: [sample1FPKM, sample2FPKM, ratio]}
    sample1 = comparisons.split(',')[0]
    sample2 = comparisons.split(',')[1]

    with open(fpkm_tracking, 'r') as fpkmfile:

    #Find out which fields you need to look at
        header = fpkmfile.readline().strip().split('\t')
        for idx, fieldname in enumerate(header):
            if fieldname == 'gene_short_name':
                genefield = idx
            if sample1 + '_FPKM' == fieldname:
                sample1FPKMfield = idx
            if sample1 + '_status' == fieldname:
                sample1statusfield = idx
            if sample2 + '_FPKM' == fieldname:
                sample2FPKMfield = idx
            if sample2 + '_status' == fieldname:
                sample2statusfield = idx

        for line in fpkmfile:
            line = line.strip().split('\t')
            if line[sample1statusfield] == 'OK' and line[sample2statusfield] == 'OK' and float(line[sample1FPKMfield]) >=fpkmcutoff and float(line[sample2FPKMfield]) >=fpkmcutoff:
                gene = line[genefield]
                sample1FPKM = float(line[sample1FPKMfield])
                sample2FPKM = float(line[sample2FPKMfield])
                #log2 of ratio of sample2 to sample1
                ratio = log(sample2FPKM / sample1FPKM, 2)
                FPKMs[gene] = [sample1FPKM, sample2FPKM, ratio]

    return FPKMs

    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fpkm_tracking', type = str, help = 'genes.fpkm_tracking from Cuffdiff.', required = True)
    parser.add_argument('--comparisons', type = str, help = 'Comma separated (no spaces) list of samples with colon between groups of samples. (e.g. N2ASoma,N2AAxon:CADSoma,CADAxon)', required = True)
    parser.add_argument('--fpkmcutoff', type = int, help = 'Only consider genes with an FPKM above this in all samples.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    args = parser.parse_args()

    if ':' in args.comparisons:
        ratioofratios = {} #{gene : [sample1AFPKM, sample1BFPKM, sample2AFPKM, sample2BFPKM, ratioA, ratioB, ratioB - ratioA]}
        group1 = args.comparisons.split(':')[0]
        group2 = args.comparisons.split(':')[1]
        group1FPKMs = getFPKMratios(args.fpkm_tracking, group1, args.fpkmcutoff)
        group2FPKMs = getFPKMratios(args.fpkm_tracking, group2, args.fpkmcutoff)
        for gene in group1FPKMs:
            if gene in group2FPKMs:
                sample1AFPKM = group1FPKMs[gene][0]
                sample1BFPKM = group1FPKMs[gene][1]
                ratioA = group1FPKMs[gene][2]
                sample2AFPKM = group2FPKMs[gene][0]
                sample2BFPKM = group2FPKMs[gene][1]
                ratioB = group2FPKMs[gene][2]
                doubleratio = ratioB - ratioA #since these are already logged, just take the difference
                ratioofratios[gene] = [sample1AFPKM, sample1BFPKM, sample2AFPKM, sample2BFPKM, ratioA, ratioB, ratioB - ratioA]

        with open(args.outfile, 'w') as outfh:
            sample1A = group1.split(',')[0]
            sample1B = group1.split(',')[1]
            sample2A = group2.split(',')[0]
            sample2B = group2.split(',')[1]
            outfh.write('gene'+'\t'+sample1A+'_FPKM'+'\t'+sample1B+'_FPKM'+'\t'+sample2A+'_FPKM'+'\t'+sample2B+'_FPKM'+'\t'+'GroupARatio'
                        +'\t'+'GroupBRatio'+'\t'+'RatioOfRatios'+'\n')
            for gene in ratioofratios:
                outfh.write(('\t').join([gene, str(ratioofratios[gene][0]), str(ratioofratios[gene][1]), str(ratioofratios[gene][2]), str(ratioofratios[gene][3]), str(ratioofratios[gene][4]), str(ratioofratios[gene][5]), str(ratioofratios[gene][6])]) + '\n')
        
    elif ':' not in args.comparisons:
        group1 = args.comparisons
        FPKMs = getFPKMratios(args.fpkm_tracking, group1, args.fpkmcutoff)
        print len(FPKMs)
