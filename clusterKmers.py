#Calculates levenshtein distances between kmers and creates matrix of those distances.  To cluster and plot dendrogram, need to bring ouptut of this script into R and do:

#library(ggdendro)
#levdistances <- read.table('output.txt', header=TRUE, sep='\t')
#ggdendrogram(hclust(dist(t(levdistances))), size = 1)

############################################
#Alternatively, here's some R code to make a prettier dendrogram
#library(ggdendro)
#library(ggplot2)
#more_loc <- read.table('output.txt', header=TRUE)
#hc <- hclust(dist(t(more_loc)))
#hcdata <- dendro_data(hc, type = 'rectangle')
#labs <- label(hcdata)
#myvector <- rep('none',nrow(labs))
#myvector[which(grepl("^[[:upper:]]+$", labs$label))] = 'Enriched'
#myvector[which(grepl("^[[:lower:]]+$", labs$label))] = 'Depleted'
#labs <- data.frame(labs, group = myvector)

#ggplot() + geom_segment(data = segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data = label(hcdata), aes(x=x, y=y, label=label, hjust = 0,
#                                                                                                                           color = labs$group), size = 3) +
#  coord_flip() + scale_y_reverse(expand = c(0.2,0)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border =
#                                                                           element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
#                                                                         axis.text.x = element_blank(), axis.text.y = element_blank()) + xlab('') + ylab('') +
#  scale_color_manual(values = c('mediumslateblue', 'tomato'), name = '') + ggtitle('Kmers in \"more localized\" distal ALEs')
##############################################
#kmers with enrichments less than 1 are output in lower case.  Those above 1 are output in upper case.

#Expects output of kmerenrichment.py as input. e.g.
#kmer	File1_count	File2_count	enrichment	fisher's exact p	BH corrected p
#AUACAG	135	121	1.84	1.148103908e-06	0.00470148550327
#UCUUCU	166	418	0.65	2.47200736951e-06	0.00506143508908
#CCUGCU	225	241	1.54	4.11535532643e-06	0.00561746002058

#Usage: python clusterKmers.py <input.txt> <bh_pvaluecutoff> <outfile>

import sys
import levenshtein
import numpy as np
#from hamming_distance import hamdist



def levenshtein_kmers(kmerfile, bh_pvalue_cutoff, outfile):
    enrichments = {}
    levenshteins = {}
    kmers = []
    
    kmersfh = open(kmerfile, 'r')

    for line in kmersfh:
        line = line.strip().split('\t')
        kmer = line[0]
        enrichment = line[3]
        bh_pvalue = line[5]
        
        #Define case by enrichment (greater than or less than 1)
        if kmer != 'kmer':
            if float(bh_pvalue) <= float(bh_pvalue_cutoff):
                if enrichment == 'NA':
                    kmers.append(kmer.upper())
                    enrichments[kmer] = enrichment
                elif float(enrichment) > 1:
                    kmers.append(kmer.upper())
                    enrichments[kmer] = enrichment
                elif float(enrichment) < 1:
                    kmers.append(kmer.lower())
                    enrichments[kmer] = enrichment
                

    kmersfh.close()

    #Compare every kmer to every other kmer.  Make a list of distances for every kmer.
    #Make sure kmers are converted to upper case so distances do not include case.
    #Put that list in a dictionary with kmer as key
    for kmer in kmers:
        levenshtein_list = []
        for KMER in kmers:
            levenshtein_list.append(levenshtein.levenshtein(kmer.upper(), KMER.upper()))

        levenshteins[kmer] = levenshtein_list

    #Turn these lists into an array
    levenshtein_array = np.array(levenshteins[kmers[0]])
    for idx,kmer in enumerate(kmers):
        if idx > 0:
            levenshtein_array = np.vstack([levenshtein_array, levenshteins[kmer]])

    print levenshtein_array
    #Output array
    np.savetxt(outfile, levenshtein_array, delimiter = '\t', fmt = '%2f')

    #Put kmer list (in order) at top of output
    with open(outfile, 'r+') as f:
        old = f.read()
        f.seek(0)
        f.write(('\t').join(kmers) + '\n' + old)

    


if __name__ == '__main__':
    levenshtein_kmers(sys.argv[1], sys.argv[2], sys.argv[3])
            
