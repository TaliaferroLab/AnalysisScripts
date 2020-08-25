import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import random

#Standalone usage: python subsampleGC.py <modelfasta> <fastatomakelooklikemodel> <outfile.txt>
#As module, returns a list of fasta records.
#Each record is itself a list where the first item is the ID and the second is the sequence

def subsampleGC(modelfasta, subsamplefasta):
    #modelbins and subsamplebins are dictionaries where key is bin (e.g. '48_to_50') and value is list of IDs
    modelbins = {}
    subsamplebins = {}
    subsampledIDs = []
    subsampledfasta = []
    modelfastarecords = 0
    subsamplefastarecords = 0

    #Populate modelbins
    for record in SeqIO.parse(modelfasta, 'fasta'):
        lowerbound = 20
        upperbound = 22
        GCcontent = float(GC(record.seq))
        modelfastarecords +=1

        while upperbound <= 70:
            if GCcontent >= lowerbound and GCcontent < upperbound:
                if '%s_to_%s' % (lowerbound, upperbound) in modelbins:
                    modelbins['%s_to_%s' % (lowerbound, upperbound)].append(record.id)
                    break
                
                else:
                    modelbins['%s_to_%s' % (lowerbound, upperbound)] = [record.id]
                    break

            else:
                lowerbound +=2
                upperbound +=2

    #Populate subsamplebins
    for record in SeqIO.parse(subsamplefasta, 'fasta'):
        lowerbound = 20
        upperbound = 22
        GCcontent = float(GC(record.seq))
        subsamplefastarecords +=1

        while upperbound <= 70:
            if GCcontent >= lowerbound and GCcontent < upperbound:
                if '%s_to_%s' % (lowerbound, upperbound) in subsamplebins:
                    subsamplebins['%s_to_%s' % (lowerbound, upperbound)].append(record.id)
                    break
                
                else:
                    subsamplebins['%s_to_%s' % (lowerbound, upperbound)] = [record.id]
                    break

            else:
                lowerbound +=2
                upperbound +=2

    #Number of records in each fasta file...used for calculating density
    modelfastarecords = float(modelfastarecords)
    subsamplefastarecords = float(subsamplefastarecords)
    
    for modelbin in modelbins:
        modelbinpop = float(len(modelbins[modelbin]))
        modelbindens = float((len(modelbins[modelbin]) / modelfastarecords))
        #Number of records to pick is density of that bin in modelfasta * number of records in subsamplefasta
        subsample_records_to_pick = int(round(modelbindens * subsamplefastarecords))

        if modelbin in subsamplebins:
            subsamplebinpop = float(len(subsamplebins[modelbin]))
            
            if subsamplebinpop > subsample_records_to_pick:
                #pick random records
                random_subsampled_IDs = random.sample(subsamplebins[modelbin], subsample_records_to_pick)
                subsampledIDs += random_subsampled_IDs
            elif subsamplebinpop <= subsample_records_to_pick:
                #pick all records
                subsampledIDs += subsamplebins[modelbin]

    #Reassemble fasta from chosen IDs
    for record in SeqIO.parse(subsamplefasta, 'fasta'):
       if record.id in subsampledIDs:
           subsampledfasta.append(['>' + str(record.id), record.seq])

    print 'There were %i records in the model fasta and %i in the fasta to be subsampled.  %i records were chosen in the sampling.' % (modelfastarecords, subsamplefastarecords, len(subsampledIDs))

    #Return a list of fasta records.  Each record is itself a list where the first item is the ID and the second is the sequence
    return subsampledfasta

    
if __name__ == '__main__':
    outfh = open(sys.argv[3], 'w')
    for entry in subsampleGC(sys.argv[1], sys.argv[2]):
        outfh.write(str(entry[0]) + '\n' + str(entry[1]) + '\n')

    outfh.close()
