from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
import os

#Usage: python filterreadlength.py <directory of fastq sequences> <trimlength>

def fastqtrimmer(trim, outfile, infile):

    trim = int(trim)
    counter = 0
    handle = open(outfile, "w")
    try:
        for title, seq, qual in FastqGeneralIterator(open(infile)) :
            counter +=1
            if counter % 1000000 == 0:
                print 'On read {0}'.format(counter)
            if len(seq) >= trim and len(seq) == len(qual):
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq[:trim], qual[:trim]))
    except ValueError:
        print 'Title and second title line don\'t match for read {0}.'.format(title)
    handle.close()

if __name__ == '__main__':
    for seqfile in os.listdir(sys.argv[1]):
        if seqfile.endswith('.fastq'):
            print 'Trimming {0}...'.format(os.path.basename(seqfile))
            outfile = os.path.basename(seqfile)[:-6] + '.trimmed.fastq'
            fastqtrimmer(sys.argv[2], outfile, seqfile)

