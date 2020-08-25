#Takes fasta file of sequences and makes histogram of GC contents

#Usage: python plotGC <sequences1.fasta> <sequences2.fasta>

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def plotmultipleLength(fasta1, fasta2):
    fasta1lengths = []
    fasta2lengths = []
    fasta3lengths = []
    fasta4lengths = []
    fasta5lengths = []
    fasta6lengths = []
    fasta7lengths = []
    fasta8lengths = []
    
    for record in SeqIO.parse(fasta1, 'fasta'):
        fasta1lengths.append(float(len(record.seq)))

    for record in SeqIO.parse(fasta2, 'fasta'):
        fasta2lengths.append(float(len(record.seq)))

    fasta1median = np.median(fasta1lengths)
    fasta2median = np.median(fasta2lengths)
    Lengths = [fasta1lengths] + [fasta2lengths]

    
    plt.hist(Lengths, 40, range = [0,4000], alpha=0.5, histtype = 'step', color = ['green', 'midnightblue'], label = ['Mito', 'Ctrl'])
    

    plt.xlabel('Length (nt)')
    plt.ylabel('Count')
    plt.title('')
    plt.legend()

    #Vertical lines for medians
    plt.plot([fasta1median, fasta1median], [0, 20], color = 'green', linestyle = '-.', linewidth = 2)
    plt.plot([fasta2median, fasta2median], [0, 20], color = 'midnightblue', linestyle = '-.', linewidth = 2)
    
    
    plt.show()
    
plotmultipleLength(sys.argv[1], sys.argv[2])
