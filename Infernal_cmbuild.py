#Will take a directory full of stockholm alignments and feed them to Infernal's cmbuild one at a time.  
#The resulting covariance models are put in sys.argv[2]

#Usage: python Infernal_cmbuild.py <directory of stockholm alignments> <directory in which to put cms>

import sys
import os
from Bio import AlignIO
import subprocess

def runcmbuild(cov_model_name, alignment):
    subprocess.check_call(['cmbuild', cov_model_name, alignment])

def make_cm_models(stockholmdirectory, cmmodeldirectory):
    stockholmdirectory = os.path.abspath(stockholmdirectory)
    alignments = [os.path.join(stockholmdirectory, alignment) for alignment in os.listdir(stockholmdirectory)]
    
    #If directory doesn't exist, make it
    if os.path.exists('./' + str(cmmodeldirectory)) == False:
        os.mkdir('./' + str(cmmodeldirectory))

    for alignment in alignments:
        cm_filename = str(os.path.basename(alignment)) + '.cm' 
        runcmbuild(cm_filename, alignment)
        #Move cm to it's proper directory
        os.rename('./' + str(cm_filename), './' + str(cmmodeldirectory) + '/' + str(cm_filename))

if __name__ == '__main__':
    make_cm_models(sys.argv[1], sys.argv[2])
