#Takes a MISO gff of ALE events and a reference gff. Reference must have stop codons explicitly annotated as such.  Searches for stop codons in each ALE.  If none, prints a warning and skips.  If >1, takes the stop codon that is most proximal to the gene.  Returns gff of 3' UTR portions of each stop-codon-containing ALE.

#ASSUMES ALE ISOFORMS ARE SINGLE EXONS.  This is true almost all of the time.

#Standalone usage: python get3UTRcoords.py <ALEIsoforms> <referencegff> <outfile>

import sys

def get3UTRcoords(ALEgff, stopcodongff):

    ALEfh = open(ALEgff, 'r')
    stopcodonfh = open(stopcodongff, 'r')
    
    stopcodons = []
    ALEs = []
    
    for line in stopcodonfh:
        line = line.strip()
        line = line.split('\t')
        if line[2] == 'stop_codon':
            stopcodons.append(line)

    for line in ALEfh:
        line = line.strip()
        line = line.split('\t')
        if line[2] == 'mRNA':
            ALEs.append(line)

    ALEfh.close()
    stopcodonfh.close()


    stoplist=[]
    stoppos = []
    UTRgffAll = []

    for ALE in ALEs:
        ALEchrm = ALE[0]
        ALEstart = int(ALE[3])
        ALEstop = int(ALE[4])
        ALEstrand = ALE[6]

        for stopcodon in stopcodons:
            stopcodonchrm = stopcodon[0]
            stopcodonstart = int(stopcodon[3])
            stopcodonstop = int(stopcodon[4])
            stopcodonstrand = stopcodon[6]
            
            #if chromosomes and strands match AND stop codon start/stop lies between ALE start/stop
            if ALEchrm == stopcodonchrm and ALEstrand == stopcodonstrand and ALEstart <= stopcodonstart and ALEstop >= stopcodonstop:
                if ALEstrand == '+': #append beginning of start codon...this will include stop codon itself, useful for checking
                    stoplist.append(stopcodonstart)
                    
                elif ALEstrand == '-': #append end of start codon, which since on the - strand, is the beginning of it
                    stoplist.append(stopcodonstop)
        #Remove duplicates
        #print len(list(set(stoplist))) 
        

        if len(stoplist) > 0:
            if ALEstrand == '+': # take most proximal stop codon to give longest UTR
                stoppos = [ALEchrm, min(stoplist)]
                UTRgff = [ALEchrm, ALE[1], '3\'UTR', str(stoppos[1] + 3), str(ALEstop), ALE[5], ALEstrand, ALE[7], ALE[8]] #THIS SHOULD NOT INCLUDE STOPCODON
                
            elif ALEstrand == '-':
                stoppos = [ALEchrm, max(stoplist)]
                UTRgff = [ALEchrm, ALE[1], '3\'UTR', str(ALEstart), str(stoppos[1] - 3), ALE[5], ALEstrand, ALE[7], ALE[8]] #THIS ALSO SHOULD NOT INCLUDE STOPCODON

            UTRgffAll.append(UTRgff)

        elif len(stoplist) == 0: # if no stop codon found
            print "WARNING: no stop codon found for " + ALE[8]

        '''#This will print chromosome, length of ALE in coding region, length of ALE in UTR
        if len(stoplist) > 0:
            if ALEstrand == '+':
                print ALEchrm + '\t' + str(stoppos[1] - ALEstart) + '\t' + str(ALEstop - stoppos[1])
            elif ALEstrand == '-':
                print ALEchrm + '\t' + str(ALEstop - stoppos[1]) + '\t' + str(stoppos[1] - ALEstart)
        '''
        stoplist = [] #when done with ALE, reset lists
        stoppos = []
        UTRgff = []

    return UTRgffAll
    

    
outfh = open(sys.argv[3], 'w')

for entry in get3UTRcoords(sys.argv[1], sys.argv[2]):
    outfh.write(('\t').join(entry) + '\n')

outfh.close()
