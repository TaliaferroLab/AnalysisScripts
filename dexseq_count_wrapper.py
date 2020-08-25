#python2

import os
import sys
import subprocess

samples = ['c9ALS_01', 'c9ALS_02', 'c9ALS_03', 'c9ALS_04', 'c9ALS_05', 'c9ALS_06', 'c9ALS_07', 'c9ALS_08', 'c9ALS_09', 'c9ALS_10', 'c9ALS_11', 'c9ALS_12', 'c9ALS_13', 'c9ALS_14', 'c9ALS_15', 'c9ALS_16', 'c9ALS_17',
'Ctrl_01', 'Ctrl_02', 'Ctrl_03', 'Ctrl_04', 'Ctrl_05', 'Ctrl_06', 'Ctrl_07', 'Ctrl_08', 'Ctrl_09', 'Ctrl_10', 'Ctrl_11', 'Ctrl_12', 'Ctrl_13', 'Ctrl_14', 'Ctrl_15', 'Ctrl_16', 'Ctrl_17', 'Ctrl_18', 'Ctrl_19', 'Ctrl_20', 'Ctrl_21', 'Ctrl_22',
'sALS_01', 'sALS_02', 'sALS_03', 'sALS_04', 'sALS_05', 'sALS_06', 'sALS_07', 'sALS_08', 'sALS_11', 'sALS_12', 'sALS_13', 'sALS_14', 'sALS_15',
'SCA3_01', 'SCA3_02', 'SCA3_03', 'SCA3_04', 'SCA3_05', 'SCA3_06', 'SCA3_07', 'SCA3_08',
'SCA36_01', 'SCA36_02', 'SCA36_03', 'SCA36_04', 'SCA36_05', 'SCA36_06']

#samples = samples[0:22]
#samples = samples[22:44]
samples = samples[44:]

def makesam(STARdir, sample):
	bam = os.path.join(STARdir, sample, sample + 'Aligned.sortedByCoord.out.bam')
	command = ['samtools', 'view', bam]
	print 'Making sam for {0}...'.format(sample)
	with open(sample + '.sam', 'w') as outfh:
		subprocess.call(command, stdout = outfh)

def runcount(sample):
	sam = sample + '.sam'
	command = ['python', '/beevol/home/taliaferro/Scripts/dexseq_count.py', '-p', 'yes', '-r', 'pos', '-s', 'no', '/beevol/home/taliaferro/data/ZackEmory/BrainSamples/HTSeq/gencodecomprehensive.v28.exoncountingbins.gff', sam, sample + 'htseqcount.txt']
	print 'Counting reads for {0}...'.format(sample)
	subprocess.call(command)
	os.remove(sam)


for sample in samples:
	makesam('/beevol/home/taliaferro/data/ZackEmory/BrainSamples/STAR', sample)
	runcount(sample)


