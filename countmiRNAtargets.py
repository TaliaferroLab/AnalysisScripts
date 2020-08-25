#Designed to find how many sites there in are in a list of genes for a particular miRNA.  Takes as input a SummaryCount table from TargetScan, a list of genes, and a particular miRNA seed sequence you want to search.  Outputs to file.

#Usage: python countmiRNAseqtargets.py Targetscanfile.txt <list of genes> <list of miRNA sequences> <outfile.txt>

import sys

def countmiRNAtargets(targetscanfile, genelist, miRNAseq):
    genecounter = 0
    targetscanlines = []
    genes = []
    genes_in_table = 0
    genes_with_sites = 0
    total_conserved_sites = 0
    conserved_8mer_sites = 0
    conserved_7mer_m8_sites = 0
    conserved_7mer_1a_sites = 0
    total_nonconserved_sites = 0
    nonconserved_8mer_sites = 0
    nonconserved_7mer_m8_sites = 0
    nonconserved_7mer_1a_sites = 0
    
    targetscanfilefh = open(targetscanfile, 'r')
    for line in targetscanfilefh:
        line = line.strip().split()
        if line[0] != 'Transcript':
            targetscanlines.append(line)
    targetscanfilefh.close()

    genelistfh = open(genelist, 'r')
    for line in genelistfh:
        line = line.strip()
        genes.append(line)
    genelistfh.close()

    for gene in genes:
        gene_in_table = False
        gene_with_site = False
        genecounter +=1
        if genecounter % 100 == 0:
            sys.stderr.write('Querying gene {0} of {1}...\n'.format(genecounter, len(genes)))
        for line in targetscanlines:
            if line[1] == gene and line[3] == '10090':
                gene_in_table = True
            if line[1] == gene and line[3] == '10090' and line[2] == miRNAseq:
                gene_with_site = True
                total_conserved_sites += int(line[4])
                conserved_8mer_sites += int(line[5])
                conserved_7mer_m8_sites += int(line[6])
                conserved_7mer_1a_sites += int(line[7])
                total_nonconserved_sites += int(line[8])
                nonconserved_8mer_sites += int(line[9])
                nonconserved_7mer_m8_sites += int(line[10])
                nonconserved_7mer_1a_sites += int(line[11])
                print gene
                break
        if gene_in_table == True:
            genes_in_table +=1
        if gene_with_site == True:
            genes_with_sites +=1
    '''
    print 'For microRNA with sequence {0}...'.format(miRNAseq)
    print 'Of {0} queried genes, found {1} in the SummaryCountsTable'.format(len(genes), genes_in_table)
    print 'Found {0} total conserved sites.'.format(total_conserved_sites)
    print 'Found {0} conserved 8mer sites.'.format(conserved_8mer_sites)
    print 'Found {0} conserved 7mer_m8 sites.'.format(conserved_7mer_m8_sites)
    print 'Found {0} conserved 7mer_1a sites.'.format(conserved_7mer_1a_sites)
    print 'Found {0} total nonconserved sites.'.format(total_nonconserved_sites)
    print 'Found {0} nonconserved 8mer sites.'.format(nonconserved_8mer_sites)
    print 'Found {0} nonconserved 7mer_m8 sites.'.format(nonconserved_7mer_m8_sites)
    print 'Found {0} nonconserved 7mer_1a sites.'.format(nonconserved_7mer_1a_sites)
    print 'Found {0} total sites in {1} genes.'.format(total_conserved_sites + total_nonconserved_sites, genes_with_sites)
    total_sites = total_conserved_sites + total_nonconserved_sites
    '''

    return genes_with_sites, total_sites, total_conserved_sites, conserved_8mer_sites, conserved_7mer_m8_sites, conserved_7mer_1a_sites, total_nonconserved_sites, nonconserved_8mer_sites, nonconserved_7mer_m8_sites, nonconserved_7mer_1a_sites

if __name__ == '__main__':
    seqlistfh = open(sys.argv[3], 'r')
    seqlist = []
    for line in seqlistfh:
        line = line.strip()
        seqlist.append(line)
    seqlistfh.close()

    outfh = open(sys.argv[4], 'w')
    outfh.close()

    for seq in seqlist:
        genes_with_sites, total_sites, total_conserved_sites, conserved_8mer_sites, conserved_7mer_m8_sites, conserved_7mer_1a_sites, total_nonconserved_sites, nonconserved_8mer_sites, nonconserved_7mer_m8_sites, nonconserved_7mer_1a_sites = countmiRNAtargets(sys.argv[1], sys.argv[2], seq)

        outfh = open(sys.argv[4], 'a')
        outfh.write(seq + '\t' + str(genes_with_sites) + '\t' + str(total_sites) + '\t' + str(total_conserved_sites) + '\t' + str(conserved_8mer_sites) + '\t' + str(conserved_7mer_m8_sites) + '\t' + str(conserved_7mer_1a_sites) + '\t' +  str(total_nonconserved_sites) + '\t' + str(nonconserved_8mer_sites) + '\t' + str(nonconserved_7mer_m8_sites) + '\t' + str(nonconserved_7mer_1a_sites) + '\n')
        outfh.close()
