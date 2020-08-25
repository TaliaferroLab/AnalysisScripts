import sys
from itertools import izip
from Bio import SeqIO

#def changebarcode(oldbarcode, newbarcode, fastq):
def changebarcode(newbarcode, fastq):
    infh = open(fastq, 'r')
    linecounter = 1
    for line in infh:
        line = line.strip()
        if linecounter == 1:
            print line
        elif linecounter == 2:
            seqbarcode = line[0:4]
            #if sum(ii == jj for ii, jj in izip(seqbarcode, oldbarcode)) >= 3:
            print newbarcode + line[4:]
            #else:
            #print line
        elif linecounter == 3:
            print line
        elif linecounter == 4:
            print line
        linecounter +=1
        if linecounter == 5:
            linecounter = 1

    infh.close()

def fixfastq(input, output):
    records = []
    for record in SeqIO.parse(input, 'fastq'):
        records.append(record)

    SeqIO.write(records, output, 'fastq')

#changebarcode(sys.argv[1], sys.argv[2], sys.argv[3])
changebarcode(sys.argv[1], sys.argv[2])
