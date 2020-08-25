import gzip
import subprocess
import sys
import pdb
import os

### sample cmd line:
#parseMISO_AFE.py /net/afterthefact/data/athma/MouseEncode/RNAseq/MISO/AFE Testis_Rep1_AFE.psi
##if submitting from inside folder:
#parseMISO_AFE.py `pwd` ${tissue}_Rep${rep}_${i}.psi

## inputs:
# (1) tissue/rep specific summary file

## outputs (bed file):
# (1) chr
# (2) AFE start
# (3) AFE end
# (4) gene;PSI;CI.low;CI.high;Assigned.Counts

#example input file: 
# /net/afterthefact/data/athma/MouseEncode/RNASeq/MISO/AFE/Testis_Rep1_AFE.psi/summary_output/summary/Testis_Rep1_AFE.psi.miso_summary')

def parseMISO_AFE(FOLDER, SUMNAME):


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
    outfile=open(FOLDER +'/'+ SUMNAME +'_TUTR.bed','w') #CHANGE THIS BASED ON EVENT TYPE

    counter=0
    while ('TRUE'):
        summaryline=summary_info.readline().split()
        if not summaryline:
            break
        counter=counter+1
        chrom=summaryline[7]
        gene=summaryline[0]
        psis=summaryline[1].split(',')
        lows=summaryline[2].split(',')
        highs=summaryline[3].split(',')
        counts=summaryline[6].split(',')
        starts=summaryline[9].split(',')
        ends=summaryline[10].split(',')
        confint = []
        for idx, high in enumerate(highs):
            confint.append(float(high) - float(lows[idx]))
        if (len(psis) == 1 and len(starts) == 2):
            psi2 =round(1-float(psis[0]),2)
            outfile.write(str(chrom)+'\t'+str(starts[0])+'\t'+str(ends[0])+'\t'+str(gene)+'\t'+str(psis[0])+'\t'+str(lows[0])+'\t'+str(highs[0])+ '\t' + str(confint[0])+'\n')
            #';'+str(counts[0])+'\n')
            #The confint for one isoform should be approximately the same as the confint for the other isoform
            outfile.write(str(chrom)+'\t'+str(starts[1])+'\t'+str(ends[1])+'\t'+str(gene)+'\t'+str(psi2)+'\t'+str('NA')+'\t'+str('NA')+ '\t' + str(confint[0])+'\n')
            #';'+str('NA')+'\n')
            if (len(psis) >= 2):
                for i in range(len(psis)):
#            count_here=counts[i].split(':')[1]
                    outfile.write(str(chrom)+'\t'+str(starts[i])+'\t'+str(ends[i])+'\t'+str(gene)+'\t'+str(psis[i])+'\t'+str(lows[i])+'\t'+str(highs[i])+ '\t' + str(confint[i])+'\n')
#';'+str(count_here)+'\n')
            
    outfile.close()

    sys.stderr.write(str(counter) + ' genes processed...' + '\n')
    sys.stderr.write('ALL DONE!')      

for directory in os.listdir(sys.argv[1]):
    if os.path.isdir(directory):
        parseMISO_AFE(os.path.abspath(sys.argv[1]), os.path.basename(directory))
    
    

