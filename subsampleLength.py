#Based on length, subsample one fasta so that it its distribution of lengths looks like another.

import sys
from Bio import SeqIO
import random
import argparse

def subsampleLength(modelfasta, subsamplefasta, minlength, maxlength, numberofbins):
    modelbins = {}
    subsamplebins = {}
    subsampledIDs = []
    subsampledfasta = []
    modelfastarecords = 0
    subsamplefastarecords = 0

    #Populate modelbins
    for record in SeqIO.parse(modelfasta, 'fasta'):
        lowerbound = minlength
        upperbound = lowerbound + ((maxlength - minlength) / numberofbins)
        seqlength = len(record.seq)
        modelfastarecords +=1

        while upperbound <= maxlength:
            if seqlength >= lowerbound and seqlength < upperbound:
                if '{0}_to_{1}'.format(lowerbound, upperbound) in modelbins:
                    modelbins['{0}_to_{1}'.format(lowerbound, upperbound)].append(record.id)
                    break
                else:
                    modelbins['{0}_to_{1}'.format(lowerbound, upperbound)] = [record.id]
                    break
            else:
                #Move bounds up by one bin
                lowerbound += (maxlength - minlength) / numberofbins
                upperbound += (maxlength - minlength) / numberofbins

    #Populate subsamplebins
    for record in SeqIO.parse(subsamplefasta, 'fasta'):
        lowerbound = minlength
        upperbound = lowerbound + ((maxlength - minlength) / numberofbins)
        seqlength = len(record.seq)
        subsamplefastarecords +=1
        while upperbound <= maxlength:
            if seqlength >= lowerbound and seqlength < upperbound:
                if '{0}_to_{1}'.format(lowerbound, upperbound) in subsamplebins:
                    subsamplebins['{0}_to_{1}'.format(lowerbound, upperbound)].append(record.id)
                    break
                else:
                    subsamplebins['{0}_to_{1}'.format(lowerbound, upperbound)] = [record.id]
                    break
            else:
                #Move bounds up by one bin
                lowerbound += (maxlength - minlength) / numberofbins
                upperbound += (maxlength - minlength) / numberofbins
        
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
                #Pick random records
                random_subsampled_IDs = random.sample(subsamplebins[modelbin], subsample_records_to_pick)
                subsampledIDs += random_subsampled_IDs
            elif subsamplebinpop <= subsample_records_to_pick:
                #Pick all records
                subsampledIDs += subsamplebins[modelbin]

    #Reassemble fasta from chosen IDs
    for record in SeqIO.parse(subsamplefasta, 'fasta'):
        if record.id in subsampledIDs:
            subsampledfasta.append(['>' + str(record.id), record.seq])

    print 'There were {0} records in the model fasta and {1} in the fasta to be subsampled. {2} records were chosen in the sampling.'.format(int(modelfastarecords), int(subsamplefastarecords), len(subsampledIDs))

    #Return a list of fasta records.
    #Each record is itself a list where the first item is the ID and the second is the sequence.
    return subsampledfasta
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--modelfasta', type = str, required=True, help = 'Fasta file to use as a model.')
    parser.add_argument('--subsamplefasta', type = str, required=True, help = 'Fasta file to subsample.')
    parser.add_argument('--minlength', type = int, required=True, help = 'For length distribution, min length to use.')
    parser.add_argument('--maxlength', type = int, required=True, help = 'For length distribution , max length to use.')
    parser.add_argument('--numberofbins', type = int, required=True, help = 'For length distribution, number of bins to use.')
    parser.add_argument('--output', type = str, required=True, help = 'Output fasta file.')
    args = parser.parse_args()

    outfh = open(args.output, 'w')
    for entry in subsampleLength(args.modelfasta, args.subsamplefasta, args.minlength, args.maxlength, args.numberofbins):
        outfh.write(str(entry[0]) + '\n' + str(entry[1]) + '\n')

    outfh.close()
