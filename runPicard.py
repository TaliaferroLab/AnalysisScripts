import subprocess
import os

samples = ['c9ALS_M01', 'c9ALS_M02', 'c9ALS_M34', 'c9ALS_M57', 'c9ALS_M06', 'c9ALS_M62', 'c9ALS_M63', 'c9ALS_M07', 'sALS_M11', 'sALS_M12', 'sALS_M13', 'sALS_M14', 'sALS_M70', 'sALS_M71', 'sALS_M72', 'sALS_M73', 'sALS_M74', 'sALS_M75', 'Ctl_M24', 'Ctl_M86', 'Ctl_M90', 'Ctl_M91', 'Ctl_M93', 'Ctl_M94', 'Ctl_M95', 'Ctl_M97']


for sample in samples:
	strand = 'SECOND_READ_TRANSCRIPTION_STRAND'
	#if sample in nostrand: #these are all nonstranded
		#strand = 'NONE'
	#if sample in secondstrand:
		#strand = 'SECOND_READ_TRANSCRIPTION_STRAND'

	bam = os.path.join('/beevol/home/taliaferro/data/ZackEmory/BrainSamples/STAR/Petrucelli', sample, '{0}Aligned.sortedByCoord.out.bam'.format(sample))
	output = os.path.join('/beevol/home/taliaferro/data/ZackEmory/BrainSamples/STAR/Petrucelli', sample, '{0}.picard.txt'.format(sample))
	refflat = '/beevol/home/taliaferro/Annotations/hg38/Gencode28/refFlat.hg38.Gencode28.txt'
	command = ['java', '-jar', '/beevol/home/taliaferro/Tools/picard.jar', 'CollectRnaSeqMetrics',
	'I={0}'.format(bam), 'O={0}'.format(output), 'REF_FLAT={0}'.format(refflat), 'STRAND={0}'.format(strand)]
	print 'Picard for {0}...'.format(sample)
	subprocess.call(command)
