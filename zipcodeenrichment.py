from Bio import SeqIO
import re
import sys


def countKmers(fasta1):
    totalseqlen = 0
    for record in SeqIO.parse(fasta1, 'fasta'):
        seq = str(record.seq.transcribe())
        totalseqlen += len(seq)

        m = re.search('GGACU(.{9,18})ACA', seq)
        if m:
            print record.id, m.groups(), len(m.groups()[0])

    print totalseqlen

if __name__ == '__main__':
    countKmers(sys.argv[1])
