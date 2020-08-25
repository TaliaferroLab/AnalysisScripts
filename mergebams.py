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

proteins = ['ENCSR815CVQ','ENCSR898NWE','ENCSR245BNJ','ENCSR081EST','ENCSR913CAE','ENCSR129RWD','ENCSR661HEL','ENCSR084SCN','ENCSR771ZZL','ENCSR942UNX','ENCSR344XID','ENCSR936VPP']
for protein in proteins:
	bam1 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/controls/{0}/STAR_v4/rep1.bam'.format(protein)
	bam2 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/controls/{0}/STAR_v4/rep2.bam'.format(protein)

	mergebam(protein, bam1, bam2)
	#mergedbam = protein + '.merged.bam'
	#indexbam(mergedbam)




