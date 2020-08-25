#Takes fasta file of sequences and makes histogram of GC contents

#Usage: python plotGC <sequences.fasta>

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def plotLength(fasta):
    Lengths = []
    for record in SeqIO.parse(fasta, 'fasta'):
        Lengths.append(float(len(record.seq)))

    print 'The average length is %f' % np.mean(Lengths)
    print 'The median length is %f' % np.median(Lengths)

    plt.hist(Lengths, 30, range = [0,3000], facecolor='green', alpha=0.75)

    plt.xlabel('Length (nt)')
    plt.ylabel('Count')
    plt.title('ALE distal UTR lengths')
    
    plt.show()

plotLength(sys.argv[1])
