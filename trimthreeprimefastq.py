from Bio import SeqIO
import sys

def trimthreeprime(fastq, number_to_trim, outfile):
    recordcounter = 0
    number_to_trim = int(number_to_trim)
    outfh = open(outfile, 'w')
    for seq_record in SeqIO.parse(fastq, 'fastq'):
        recordcounter +=1
        if recordcounter % 1000000 == 0:
            sys.stderr.write('Trimming read {0}.\n'.format(recordcounter))
            
        SeqIO.write(seq_record[:-number_to_trim], outfh, 'fastq')
    sys.stderr.write('Done!')
    outfh.close()

if __name__ == '__main__':
    trimthreeprime(sys.argv[1], sys.argv[2], sys.argv[3])
