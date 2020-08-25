#Takes a .miso output file and flips the psi values for the two isoforms.  Requires that each .miso file has only two possible isoforms.  If this is not the case, the file is skipped. If it is, the comma separated values under sampled_psi are flip and the file is given a .miso.psiflip name.

#Usage: python <input dir>

#Output will be in './PSIFlip/inputdir'. 

from itertools import islice
import os
import sys

indir = sys.argv[1]

if not os.path.exists('./PSIFlip/' + indir): #If ../PSIFlip directory doesn't exist, create it
    os.makedirs('./PSIFlip/' + indir)
    

for infile in os.listdir(indir):
    if not infile.endswith('miso'): continue
    if infile.endswith('miso.psiflip'): continue
    
    outfile = infile.replace('miso', 'miso.psiflip')

    
    print 'working on', infile

    infh = open(os.path.join(indir, infile),'r')
    
    head = list(islice(infh,2)) #get first two lines of file (header) in list form

    psivalues = []

    for line in infh: #islice removed header so every line now is psivalues
        line = line.strip()
        line = line.replace(',','\t')
        line = line.split('\t')
        psivalues.append(line)

        
    if len(psivalues[0]) == 3:
        outfh = open(os.path.join('./PSIFlip/' + indir, outfile), 'w') #if two isoforms do these things
        for line in head:
            outfh.write(line)
        for psivalue in psivalues:
            outfh.write(psivalue[1] + ',' + psivalue[0] + '\t' + psivalue[2] + '\n')

        outfh.close()
    infh.close()

