#Will take directory full of calibrated covariance models and feed them to cmsearch one at a time.
#The resulting outputs are put in sys.argv[3]

#Usage: python Infernal_cmsearch.py <directory of calibrated models> <fasta file to search through> <output dir>

import sys
import os
import subprocess
        
def runcmsearch(calibratedmodel, fastafile):
    sys.stderr.write('Searching for instances of {0} in {1}.\n'.format(str(calibratedmodel), str(fastafile)))
    #Search top (+) strand only
    #Stdout is returned as a string
    cmsearchoutput = subprocess.Popen(['cmsearch', '--toponly', calibratedmodel, fastafile], stdout=subprocess.PIPE).communicate()[0]
    return cmsearchoutput

def cmsearchondirectory(calibratedmodelsdirectory, fastafile, outdirectory):
    calibratedmodelsdirectory = os.path.abspath(calibratedmodelsdirectory)
    #Models are ONLY those items in this directory that are FILES, not other things (other directories, etc.)
    calibratedmodels = [os.path.join(calibratedmodelsdirectory, calibratedmodel) for calibratedmodel in os.listdir(calibratedmodelsdirectory) if os.path.isfile(os.path.join(calibratedmodelsdirectory, calibratedmodel))]

    #if outdirectory doesn't exist, make it
    if os.path.exists('./' + str(outdirectory)) == False:
        os.mkdir('./' + str(outdirectory))

    for calibratedmodel in calibratedmodels:
        search_filename = str(os.path.basename(calibratedmodel)).replace('aln.c.cm', 'csearch')
        cmsearchoutput = runcmsearch(calibratedmodel, fastafile)
        outfh = open(search_filename, 'w')
        outfh.write(cmsearchoutput)
        outfh.close()

        #Move csearch file to outdirectory
        os.rename('./' + str(search_filename), './' + str(outdirectory) + '/' + str(search_filename))
    
if __name__ == '__main__':
    cmsearchondirectory(sys.argv[1], sys.argv[2], sys.argv[3])
