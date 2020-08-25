from Bio import SeqIO
import sys

outputhandle1 = open(sys.argv[2],'w')
outputhandle2 = open(sys.argv[3],'w')
outputhandle3 = open(sys.argv[4],'w')
outputhandle4 = open(sys.argv[5],'w')
outputhandle_nomatch = open(sys.argv[6],'w')

barcode1count = 0
barcode2count = 0
barcode3count = 0
barcode4count = 0
recordcount = 0

#allows one mismatch, can allow any number of mismatches
#barcodes at positions -8 to -3 of first line of fastq record
#Barcode 1 = GGCATG, Barcode 2 = TTATAG, Barcode 3 = AATCGT, Barcode 4 = CCGGCT
#Barcodes in fastq file are rev comp of those above
#usage: python BarcodeSplitter.py <inputfile> <Barcode1Output> <Barcode2Output> <Barcode3Output> <Barcode4Output> <NoMatchOutput>

#add more barcodes if necessary

for record in SeqIO.parse(open(sys.argv[1], 'rU'),'fastq'):
    if len(record.seq) < 29:  #only consider sequences of at least a particular length
        continue

    else:
        if (record.id[-8] == 'C'):  #first barcode
            barcode1count +=1
        if (record.id[-7] == 'A'):
            barcode1count +=1    
        if (record.id[-6] == 'T'):
            barcode1count +=1
        if (record.id[-5] == 'G'):
            barcode1count +=1
        if (record.id[-4] == 'C'):  
            barcode1count +=1
        if (record.id[-3] == 'C'):  
            barcode1count +=1    
        if (record.id[-8] == 'C'):  #second barcode
            barcode2count +=1
        if (record.id[-7] == 'T'):
            barcode2count +=1    
        if (record.id[-6] == 'A'):
            barcode2count +=1
        if (record.id[-5] == 'T'):
            barcode2count +=1
        if (record.id[-4] == 'A'):  
            barcode2count +=1
        if (record.id[-3] == 'A'):  
            barcode2count +=1
        if (record.id[-8] == 'A'):  #third barcode
            barcode3count +=1
        if (record.id[-7] == 'C'):
            barcode3count +=1
        if (record.id[-6] == 'G'):  
            barcode3count +=1
        if (record.id[-5] == 'A'):  
            barcode3count +=1
        if (record.id[-4] == 'T'):
            barcode3count +=1
        if (record.id[-3] == 'T'):
            barcode3count +=1
        if (record.id[-8] == 'A'):  #fourth barcode
            barcode4count +=1
        if (record.id[-7] == 'G'):  
            barcode4count +=1
        if (record.id[-6] == 'C'):  
            barcode4count +=1
        if (record.id[-5] == 'C'):
            barcode4count +=1    
        if (record.id[-4] == 'G'):
            barcode4count +=1
        if (record.id[-3] == 'G'):
            barcode4count +=1

    if barcode1count > 4:
        SeqIO.write(record[:], outputhandle1, 'fastq')  #if you want to remove the barcode, do it here by slicing
    elif barcode2count > 4:
        SeqIO.write(record[:], outputhandle2, 'fastq')
    elif barcode3count > 4:
        SeqIO.write(record[:], outputhandle3, 'fastq')
    elif barcode4count > 4:
        SeqIO.write(record[:], outputhandle4, 'fastq')
    else:
        SeqIO.write(record[:], outputhandle_nomatch, 'fastq')

    barcode1count = 0
    barcode2count = 0
    barcode3count = 0
    barcode4count = 0
    recordcount +=1
    if recordcount % 1000000 == 0:
        print 'Processing read %s' % recordcount

outputhandle1.close()
outputhandle2.close()
outputhandle3.close()
outputhandle4.close()
outputhandle_nomatch.close()
