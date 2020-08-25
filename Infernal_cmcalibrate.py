#Will take a directory full of covariance models and feed them to Infernal's cmcalibrate one at a time.
#The resulting covariance models are put in sys.argv[3].
#This script may take several days to complete.  Calibrating models of large alignments is slow.

#Necessary modules: biopython, infernal

#Usage: python Infernal_cmcalibrate.py <directory of cov models> <directory of stockholm alignments> <directory in which to put calibrated models>

import sys
import os
from Bio import AlignIO
import subprocess
import shutil

def runcmcalibrate(modelname):
    subprocess.check_call(['cmcalibrate', modelname])

def calibratemodels(modelsdirectory, stockholmdirectory, calibratedmodelsdirectory):
    modelcounter = 0
    calibratedmodelcounter = 0
    modelsdirectory = os.path.abspath(modelsdirectory)
    stockholmdirectory = os.path.abspath(stockholmdirectory)
    models = [os.path.join(modelsdirectory, model) for model in os.listdir(modelsdirectory)]

    #If calibratedmodelsdirectory doesn't exist, make it
    if os.path.exists ('./' + str(calibratedmodelsdirectory)) == False:
        os.mkdir('./' + str(calibratedmodelsdirectory))

    for model in models:
        modelcounter +=1
        if modelcounter % 10 == 0:
            sys.stderr.write('Calibrating model {0} of {1}.\n'.format(modelcounter, len(models)))
        #Get corresponding stockholm alignment
        #if calibrated model doesn't already exist in calibratedmodelsdirectory
        if os.path.exists(str(os.path.abspath(calibratedmodelsdirectory) + '/' + 
                              str(os.path.basename(model.replace('.cm', '.c.cm'))))) == False:
            stockholmalignmentname = os.path.basename(model.replace('.cm',''))
            stockholmalignment = os.path.join(stockholmdirectory, stockholmalignmentname)

            #Only calibrate models where the alignment is less than 2 kb.
            #Otherwise it takes forEVER. Still, expect an average of hours per alignment.
            alignment_length = AlignIO.read(stockholmalignment, 'stockholm').get_alignment_length()
            if alignment_length >= 10 and alignment_length <= 4000:
                runcmcalibrate(model)
                #Move calibrated model to its correct directory
                shutil.copy2(str(model), str(os.path.abspath(calibratedmodelsdirectory)) + '/' + 
                             str(os.path.basename(model.replace('.cm','.c.cm'))))
                calibratedmodelcounter +=1

    sys.stderr.write('Calibrated {0} of {1} models. The rest were too big.\n'.format(calibratedmodelcounter, len(models)))
        

if __name__ == '__main__':
    calibratemodels(sys.argv[1], sys.argv[2], sys.argv[3])
