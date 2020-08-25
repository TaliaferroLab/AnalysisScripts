import gzip
import subprocess
import sys
import pdb
import re
import numpy

## sample cmd line:
#parseMISO_SE.py /net/utr/data/atf/athma/MouseEncode/GingerasLongRead/MISO/SE Testis_Rep1_SE.psi
## if submitting from inside folder:
#parseMISO_SE.py `pwd` ${tissue}_Rep${rep}_${i}.psi

## inputs:
# (1) tissue/rep specific summary file

## outputs (bed file):
# (1) chr
# (2) SE start
# (3) SE end
# (4) gene;PSI;CI.low;CI.high;Assigned.Counts

# example input file:
# /net/utr/data/atf/athma/MouseEncode/GingerasLongRead/MISO/SE/Testis_Rep1_SE.psi/summary_output/summary/Testis_Rep1_AFE.psi.miso_summary

FOLDER=sys.argv[1]
# FOLDER=('/net/utr/data/atf/athma/MouseEncode/GingerasLongRead/MISO/SE')
SUMNAME=sys.argv[2]
# SUMNAME=('Testis_Rep1_SE.psi')

sys.stderr.write('This is the error for '+ FOLDER +'/'+SUMNAME+'\n')

# Open folder and flush header line
summary_info=open(FOLDER+'/'+SUMNAME+'/summary_output/summary/'+SUMNAME+'.miso_summary')
summaryinfo=summary_info.readline().split()

##### Example lines
###['event_name', 'miso_posterior_mean', 'ci_low', 'ci_high', 'isoforms', 'counts', 'assigned_counts', 'chrom', 'strand', 'mRNA_starts', 'mRNA_ends']
## event name
#['13:59710778:59710954:+:ENSMUSG00000021555',
## miso_posterior_mean
#'0.02,0.01,0.03,0.01,0.90,0.02',
## 'ci_low'
#'0.00,0.00,0.00,0.00,0.83,0.00',
## 'ci_high'
#'0.05,0.04,0.07,0.07,0.95,0.05',
##isoforms
#'13:59710778:59710954:+:ENSMUSG00000021555.1.A.up_13:59710778:59710954:+:ENSMUSG00000021555.1.A.sk_13:59710778:59710954:+:ENSMUSG00000021555.1.A.dn',
#'13:59710778:59710954:+:ENSMUSG00000021555.1.B.up_13:59710778:59710954:+:ENSMUSG00000021555.1.B.dn',
#'13:59710778:59710954:+:ENSMUSG00000021555.2.A.up_13:59710778:59710954:+:ENSMUSG00000021555.2.A.sk_13:59710778:59710954:+:ENSMUSG00000021555.2.A.dn',
#'13:59710778:59710954:+:ENSMUSG00000021555.2.B.up_13:59710778:59710954:+:ENSMUSG00000021555.2.B.dn',
#'13:59710778:59710954:+:ENSMUSG00000021555.3.A.up_13:59710778:59710954:+:ENSMUSG00000021555.3.A.sk_13:59710778:59710954:+:ENSMUSG00000021555.3.A.dn',
#'13:59710778:59710954:+:ENSMUSG00000021555.4.A.up_13:59710778:59710954:+:ENSMUSG00000021555.4.A.sk_13:59710778:59710954:+:ENSMUSG00000021555.4.A.dn'
##counts
#'(0,0,0,0,0,0):209,(0,0,0,0,1,0):131,(1,0,0,0,1,0):7,(1,0,1,0,0,0):8,(1,0,1,0,1,1):21',
## assigned counts
#'0:0,1:0,2:9,3:0,4:158',
## chrom
#'chr13',
## strand
#'+',
## mRNA starts
#'59710122,59710122,59709506,59709506,59710122,59709506',
## mRNA ends
#'59714176,59714176,59714176,59714176,59714176,59714176'

# Open output folder
outfile=open(FOLDER +'/'+ SUMNAME +'.bed','w')

def sktxpts(txlist):
    txhold=[s for s in txlist if "sk" in s]
    hold=[]
    for x in txhold:
        hold.append(txlist.index(x))
    return(hold)
        
counter=0
while('TRUE'):
    summaryline=summary_info.readline().split()
    if not summaryline:
        break
    counter=counter+1
    chrom=summaryline[7]
    gene=summaryline[0].split(':')[4]
    start=summaryline[0].split(':')[1]
    end=summaryline[0].split(':')[2]
    psis=numpy.array(summaryline[1].split(','),dtype='float')
    counts=numpy.array([f.split('):')[1] for f in summaryline[5].split(',(')][1:],dtype='int').sum()
    txpts=summaryline[4].split(',')
    sk_psi=round(psis[sktxpts(txpts)].sum(),3)
    outfile.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(gene)+';'+str(sk_psi)+';'+str(counts)+'\n')

outfile.close()

sys.stderr.write(str(counter) + ' genes processed...' + '\n')
sys.stderr.write('ALL DONE!')
    
    
