#Takes a .miso_bf file from run_miso.py --compare-samples.  
#These files contain events with more than two isoforms. 
#This script removes those events.

#Usage: python MISO_2Isoforms.py <infile> <outfile>

import sys


def getTwoIsoformEvents(miso_bffile, outfile):
    twoisoformevents = []
    infh = open(miso_bffile, 'r')
    for line in infh:
        line = line.strip().split('\t')
        #This is the header line
        if line[0] == 'event_name':
            header = line
        #Multiple isoform events have psi values separated by commas
        elif ',' not in line[1]:
            twoisoformevents.append(line)

    infh.close()
    outfh = open(outfile, 'w')
    outfh.write(('\t').join(header) + '\n')
    for event in twoisoformevents:
        outfh.write(('\t').join(event) + '\n')

    outfh.close()
    
if __name__ == '__main__':
    getTwoIsoformEvents(sys.argv[1], sys.argv[2])
