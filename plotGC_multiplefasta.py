#Takes fasta file of sequences and makes histogram of GC contents

#Usage: python plotGC_multiplefasta.py <sequences1.fasta> <sequences2.fasta>

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def plotGC_multiplefasta(fasta1, fasta2):
    fasta1GCs = []
    fasta2GCs = []
    
    
    
    for record in SeqIO.parse(fasta1, 'fasta'):
        fasta1GCs.append(float(GC(record.seq[:]))) #slice

    for record in SeqIO.parse(fasta2, 'fasta'):
        fasta2GCs.append(float(GC(record.seq[:]))) #slice

    fasta1median = np.median(fasta1GCs)
    fasta2median = np.median(fasta2GCs)
    
    GCs = [fasta1GCs] + [fasta2GCs]

    plt.hist(GCs, 25, range = [20,70], alpha=0.5, histtype = 'bar', color = ['blue', 'green'], label=['fasta1','fasta2'])

    plt.xlabel('GC percentage')
    plt.ylabel('Count')
    plt.title('GC content of UTRs excluding last 50 nt')
    plt.legend()

    #Vertical lines for medians
    plt.plot([fasta1median, fasta1median], [0, 60], color = 'blue', linestyle = '-.', linewidth = 2)
    plt.plot([fasta2median, fasta2median], [0, 60], color = 'green', linestyle = '-.', linewidth = 2)
    
    
    
    plt.show()

plotGC_multiplefasta(sys.argv[1], sys.argv[2])
