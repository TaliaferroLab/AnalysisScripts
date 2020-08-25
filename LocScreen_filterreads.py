from Bio import SeqIO
import sys
import os

def filterreads(reads):
    filteredreads = []
    analyzedreads = 0
    output = os.path.basename(reads).split('.')[0] + '.filtered.fasta'
    
    for record in SeqIO.parse(reads, 'fastq'):
        analyzedreads +=1
        if analyzedreads % 1000000 == 0:
            print 'Analyzing read {0}...'.format(analyzedreads)
        if record.seq.endswith('TGGAATTCTCGGGTGCCAAGG'):
            filteredreads.append(record)

    print 'Analyzed {0} reads and kept {1} of them.'.format(analyzedreads, len(filteredreads))

    SeqIO.write(filteredreads, output, 'fasta')

filterreads(sys.argv[1])
