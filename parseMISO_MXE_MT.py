import gzip
import subprocess
import sys
import pdb
import os

### sample cmd line:
#parseMISO_MXE.py /net/afterthefact/data/athma/MouseEncode/RNAseq/MISO/AFE Testis_Rep1_AFE.psi
##if submitting from inside folder:
#parseMISO_MXE.py `pwd` ${tissue}_Rep${rep}_${i}.psi

## inputs:
# (1) tissue/rep specific summary file

## outputs (bed file):
# (1) chr
# (2) AFE start
# (3) AFE end
# (4) gene;PSI;CI.low;CI.high;Assigned.Counts

#example input file: 
# /net/afterthefact/data/athma/MouseEncode/RNASeq/MISO/AFE/Testis_Rep1_AFE.psi/summary_output/summary/Testis_Rep1_AFE.psi.miso_summary')

def parseMISO_MXE(FOLDER, SUMNAME):


    #FOLDER=('/net/afterthefact/data/athma/MouseEncode/RNAseq/MISO/AFE')
    
    #SUMNAME=('Testis_Rep1_AFE.psi')
    #SUMNAME=('Testis_Rep1_AFE.psi/summary_output/summary/Testis_Rep1_AFE.psi.miso_summary')

    sys.stderr.write('This is the error for ' + FOLDER +'/'+ SUMNAME +'\n')
    
    # Open folder and flush header line
    summary_info=open(FOLDER+'/'+SUMNAME+'/summary_output/summary/'+SUMNAME+'.miso_summary')
    summaryinfo=summary_info.readline().split()

    ### Examples lines
    #ENSMUSG00000066151@4:61961182:62021610	0.70	0.24	0.99	'ENSMUSG00000066151@4:61961182:62021610.A.0','ENSMUSG00000066151@4:61961182:62021610.B.0'	(0,0):28,(1,0):1	0:1	chr4	-	62021404,61965225	62021610,61966057

    # Open output folder
    outfile=open(FOLDER +'/'+ SUMNAME +'_MXE.bed','w') #CHANGE THIS BASED ON EVENT TYPE

    counter=0
    while ('TRUE'):
        summaryline=summary_info.readline().split()
        if not summaryline:
            break
        counter=counter+1
        chrom=summaryline[7]
        gene=summaryline[0]
        #4 psi values here.  
        psis=summaryline[1].split(',')
        lows=summaryline[2].split(',')
        highs=summaryline[3].split(',')
        counts=summaryline[6].split(',')
        starts=summaryline[9].split(',')
        ends=summaryline[10].split(',')
        #Only consider lines with 4 psi values.  Only having one psi value doesn't make sense.
        if len(psis) == 4:
            confint = []
            for idx, high in enumerate(highs):
                confint.append(float(high) - float(lows[idx]))
                
            outfile.write(('\t').join([chrom, starts[0], ends[0], gene + '_none', psis[0], lows[0], highs[0], str(confint[0])]) + '\n')
            outfile.write(('\t').join([chrom, starts[1], ends[1], gene + '_mxe1', psis[1], lows[1], highs[1], str(confint[1])]) + '\n')
            outfile.write(('\t').join([chrom, starts[2], ends[2], gene + '_mxe2', psis[2], lows[2], highs[2], str(confint[2])]) + '\n')
            outfile.write(('\t').join([chrom, starts[3], ends[3], gene + '_both', psis[3], lows[3], highs[3], str(confint[3])]) + '\n')
            
            
    outfile.close()

    sys.stderr.write(str(counter) + ' genes processed...' + '\n')
    sys.stderr.write('ALL DONE!')      

for directory in os.listdir(sys.argv[1]):
    if os.path.isdir(directory):
        parseMISO_MXE(os.path.abspath(sys.argv[1]), os.path.basename(directory))
    
    

