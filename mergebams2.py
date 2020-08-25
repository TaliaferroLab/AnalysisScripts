import subprocess
import os
from datetime import datetime

def mergebam(protein, bam1, bam2):
	print 'Merging bams for {0}...'.format(protein)
	print str(datetime.now())
	subprocess.call(['samtools', 'merge', protein + '.merged.bam', bam1, bam2])
	print 'Done merging bams for {0}!'.format(protein)

def indexbam(bam):
	print 'Indexing {0}...'.format(bam)
	print str(datetime.now())
	subprocess.call(['samtools', 'index', bam])
	print 'Done indexing {0}!'.format(bam)

proteins = ['NONO','PA2G4','PAPOLA','PCBP2','PES1','POLR2G','PPIL4','PRPF6','PRPF8','RAVER1','RBM15','RBM22','RBM3','RBM34','RPS10','RPS19','RPS3','RRP9','RTF1','SERBP1','SF1','SLTM','SMNDC1','SND1','SRP68','SRSF4','SRSF5','SRSF7','SRSF9','SSB','SUB1','SUPV3L1','TAF15','TARDBP','TBRG4','TFIP11','TIA1','TRA2A','U2AF1','U2AF2','UCHL5','XRCC5','XRCC6','XRN2','ZRANB2']
for protein in proteins:
	bam1 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/{0}/STAR_v4/rep1.bam'.format(protein)
	bam2 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/{0}/STAR_v4/rep2.bam'.format(protein)

	mergebam(protein, bam1, bam2)
	#mergedbam = protein + '.merged.bam'
	#indexbam(mergedbam)