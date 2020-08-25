#Trim all reads longer than <length> to that length.  If a read is shorter, it still passes and is written.

from Bio import SeqIO
import argparse
import sys

def trimreads(fastq, length, fastqout):
    length = int(length)
    recordcounter = 0
    for record in SeqIO.parse(fastq, 'fastq'):
        recordcounter +=1
        if recordcounter % 10000000 == 0:
            sys.stderr.write('Considering read {0}.\n'.format(recordcounter))
        if len(record.seq) >= length:
            with open(fastqout, 'a') as fastqoutfh:
                SeqIO.write(record[:length], fastqoutfh, 'fastq')
        elif len(record.seq) < length:
            with open(fastqout, 'a') as fastqoutfh:
                SeqIO.write(record, fastqoutfh, 'fastq')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq', type = str, help = 'Fastq of reads.', required = True)
    parser.add_argument('--maxlength', type = str, help = 'Max length of output reads.', required = True)
    parser.add_argument('--fastqout', type = str, help = 'Outfile for trimmed reads.', required = True)
    args = parser.parse_args()

    fastqoutfh = open(args.fastqout, 'w')
    fastqoutfh.close()
    trimreads(args.fastq, args.maxlength, args.fastqout)
