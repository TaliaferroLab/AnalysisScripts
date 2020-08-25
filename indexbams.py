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

proteins = ['AATF','AGO3','BCCIP','BCLAF1','CELF1','CNOT8','CPSF6','CSTF2','CSTF2T','DAZAP1','DDX21','DDX24','DDX27','DDX28','DDX47','DDX51','DDX52','DDX55','EEF2','EFTUD2','EIF2C1','EIF3D','EIF3G','EWSR1','EXOSC9','FAM120A','FASTKD2','FUBP3','FXR2','GEMIN5','GTF2F1','HLTF','HNRNPC','HNRNPU','IGF2BP1','NFX1']
for protein in proteins:
	bam1 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/{0}/STAR_v4/rep1.bam'.format(protein)
	bam2 = '/net/uorf/data/backup/pfreese/encode_KD/raw_data/k562/{0}/STAR_v4/rep2.bam'.format(protein)

	#mergebam(protein, bam1, bam2)
	mergedbam = protein + '.merged.bam'
	if os.path.isfile('/net/nevermind/data/nm/jmtali/ENCODEkds/STARruns/K562/bams/{0}.merged.bam.bai'.format(protein)):
		print 'Index already exists for {0}!'.format(protein)
		continue
	else:
		indexbam(mergedbam)




