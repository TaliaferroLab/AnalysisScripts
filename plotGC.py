#Takes fasta file of sequences and makes histogram of GC contents

#Usage: python plotGC <sequences.fasta>

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def plotGC(fasta):
    GCs = []
    for record in SeqIO.parse(fasta, 'fasta'):
        GCs.append(float(GC(record.seq[:]))) #slice

    print np.mean(GCs)
    print np.median(GCs)

    plt.hist(GCs, 25, range = [20,70], facecolor='green', alpha=0.75)

    plt.xlabel('GC percentage')
    plt.ylabel('Count')
    plt.title('ALE distal UTR %GC')

    plt.show()

plotGC(sys.argv[1])
