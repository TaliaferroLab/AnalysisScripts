#Takes a gff of the genome, a list of genes you are interested in, and fasta files of the genome.  Returns amino acid sequences for the longest ORF of each gene of interest. Genes of interest file has Ensembl gene names in second field of tab-separated file

#Usage: python longestORF.py <gff> <genes of interest> <fasta files directory> <outfile>

import sys
import os
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from operator import itemgetter

def longestORF(gff, genes_of_interest_file, sequencefile):
    longestORFs = {} #{[geneID, strand] : [[exon1start, exon1stop], ... , [exonNstart, exonNstop]] for longest CDS in gene}
    genes_of_interest = []
    seqs = {} #dictionary where key is ID and value is sequence
    genes_with_ORFs = 0
    total_ORFs = 0
    counter = 0
    
    #Get genes of interest
    genesfh = open(genes_of_interest_file, 'r')
    for line in genesfh:
        line = line.strip().split('\t')
        genes_of_interest.append(line[1])
    genesfh.close()
    
    #Make gff database
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)

    genes = db.features_of_type('gene')

    #Make sequence database
    sys.stderr.write('Creating sequence dictionary...\n')
    seq_dict = SeqIO.to_dict(SeqIO.parse(sequencefile, 'fasta'))

    #find longest ORFs
    for gene in genes:
        ORFlengths = {}
        longestCDS = []
        for isoform in db.children(gene, featuretype = 'mRNA'):
            if sum(len(i) for i in db.children(isoform, featuretype = 'CDS')) > 0:
                ORFlengths[isoform.id] = sum(len(i) for i in db.children(isoform, featuretype = 'CDS'))
        if len(ORFlengths) > 0:
            ORFlengths = sorted(ORFlengths.items(), key = lambda x: x[1], reverse=True) #now a list of tuples
            longest_CDS_transcript = ORFlengths[0][0]
            total_ORFs += len(ORFlengths)
            genes_with_ORFs +=1

            for isoform in db.children(gene, featuretype = 'mRNA'):
                if isoform.id == longest_CDS_transcript:
                    for CDS in db.children(isoform, featuretype = 'CDS'):
                        longestCDS.append([CDS.start, CDS.stop])
                        longestCDS = sorted(longestCDS, key = itemgetter(0)) #sort CDS exons by start pos
                
        longestORFs[gene.id] = longestCDS

    sys.stderr.write('Found {0} total ORFs. {1} genes had at least one ORF. \n'.format(total_ORFs, genes_with_ORFs))
    genes = db.features_of_type('gene') #remake generator

    for gene in genes:
        if gene.id in genes_of_interest:
            counter +=1
            sequence = ''
            if gene.strand == '+':
                for exon in longestORFs[gene.id]:
                    exonstart = int(exon[0])
                    exonstop = int(exon[1])
                    sequence += seq_dict[gene.chrom].seq[exonstart-1:exonstop]

            elif gene.strand == '-':
                for exon in reversed(longestORFs[gene.id]): #reverse exon order
                    exonstart = int(exon[0])
                    exonstop = int(exon[1])
                    sequence += seq_dict[gene.chrom].seq[exonstart-1:exonstop].reverse_complement()

            ORFsequence = Seq(str(sequence), IUPAC.unambiguous_dna)
                
            seqs[gene.id] = ORFsequence.translate()

            if counter <= 50 and counter % 10 == 0:
                sys.stderr.write('Retrieving sequence %i \n' % (counter))
            elif counter > 50 and counter % 50 == 0:
                sys.stderr.write('Retrieving sequence %i \n' % (counter))

    sys.stderr.write('Retrieved {0} sequences.'.format(len(seqs)))

    '''
    os.remove(db_fn)
    '''

    return seqs


if __name__ == '__main__':
    seqdictionary = longestORF(sys.argv[1], sys.argv[2], sys.argv[3])
    outfh = open(sys.argv[4], 'w')
    for ID in seqdictionary:
        outfh.write('>' + ID + '\n' + str(seqdictionary[ID]) + '\n')

    outfh.close()
