from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
import argparse

def addAdapters(fasta, adapter5, adapter3, outfile):
    adapter5 = adapter5.upper()
    adapter3 = adapter3.upper()
    outseqs = []
    for record in SeqIO.parse(fasta, 'fasta'):
        newseq = adapter5 + str(record.seq).strip() + adapter3
        with open(outfile, 'a') as f:
            f.write('>' + record.id + '\n' + newseq + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type = str, help = 'Input fasta.')
    parser.add_argument('--adapter5', type = str, help = 'Adapter sequence to add to 5\' end.')
    parser.add_argument('--adapter3', type = str, help = 'Adapter sequence to add to 3\' end.')
    parser.add_argument('--output', type = str, help = 'Output file.')
    args = parser.parse_args()

    addAdapters(args.fasta, args.adapter5, args.adapter3, args.output)
