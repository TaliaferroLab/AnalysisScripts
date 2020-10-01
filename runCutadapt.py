#python3
import subprocess

samples = ['CAD_Neurite_FF_1_S13_L004_R1_001', 'CAD_Neurite_FF_2_S14_L004_R1_001', 'CAD_Neurite_FF_3_S15_L004_R1_001', 'CAD_Neurite_FF_4_S16_L004_R1_001',
'CAD_Soma_FF_1_S9_L004_R1_001', 'CAD_Soma_FF_2_S10_L004_R1_001', 'CAD_Soma_FF_3_S11_L004_R1_001', 'CAD_Soma_FF_4_S12_L004_R1_001',
'N2A_Neurite_FF_1_S5_L004_R1_001', 'N2A_Neurite_FF_2_S6_L004_R1_001', 'N2A_Neurite_FF_3_S7_L004_R1_001', 'N2A_Neurite_FF_4_S8_L004_R1_001',
'N2A_Soma_FF_1_S1_L004_R1_001', 'N2A_Soma_FF_2_S2_L004_R1_001', 'N2A_Soma_FF_3_S3_L004_R1_001', 'N2A_Soma_FF_4_S4_L004_R1_001',
'CAD_Neurite_GFP_1_S29_L004_R1_001', 'CAD_Neurite_GFP_2_S30_L004_R1_001', 'CAD_Neurite_GFP_3_S31_L004_R1_001', 'CAD_Neurite_GFP_4_S32_L004_R1_001',
'CAD_Soma_GFP_1_S25_L004_R1_001', 'CAD_Soma_GFP_2_S26_L004_R1_001', 'CAD_Soma_GFP_3_S27_L004_R1_001', 'CAD_Soma_GFP_4_S28_L004_R1_001',
'N2A_Neurite_GFP_1_S21_L004_R1_001', 'N2A_Neurite_GFP_2_S22_L004_R1_001', 'N2A_Neurite_GFP_3_S23_L004_R1_001', 'N2A_Neurite_GFP_4_S24_L004_R1_001',
'N2A_Soma_GFP_1_S17_L004_R1_001', 'N2A_Soma_GFP_2_S18_L004_R1_001', 'N2A_Soma_GFP_3_S19_L004_R1_001', 'N2A_Soma_GFP_4_S20_L004_R1_001']


#adapter1 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter2 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'
adapter1 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
adapter2 = 'CTGATCAGCGGGTTTCACTAGTGCGACCGCAAGAG'
adapter3 = 'TGGTGGCTGGTGTGGCCAAGCTTCGATATCCGCATGCTA'
adapter4 = 'GGGAAAAAGATCTCAGTGGTATTTGTGAGCCAGCACTAGTGCGACCGCAAGAG'
#adapter5 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter6 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'


for idx, sample in enumerate(samples):
	print('Trimming {0}, sample {1} of {2}...'.format(sample, idx + 1, len(samples)))
	input1 = sample + '.fastq.gz'
	input2 = input1.replace('_R1_', '_R2_')
	if 'Soma_FF' in input1:
		index = 13
	elif 'Neurite_FF' in input1:
		index = 16
	elif 'Soma_GFP' in input1:
		index = 14
	elif 'Neurite_GFP' in input1:
		index = 17
	
	samplename = sample[:index]
	samplename = samplename[:-2] + '_Rep' + samplename[-1]
	output1 = '{0}.R1.trimmed.fq.gz'.format(samplename)
	output2 = '{0}.R2.trimmed.fq.gz'.format(samplename)
	statsout = '{0}.cutadaptstats.txt'.format(samplename)
	if 'FF' in samplename:
		adapters = [adapter1, adapter2]
	elif 'GFP' in samplename:
		adapters = [adapter3, adapter4]
	command = ['cutadapt', '-g', adapters[0], '-G', adapters[1], '--minimum-length', '75', '-j', '8',
	'-o', output1, '-p', output2, input1, input2]
	#command = ['cutadapt', '-a', adapter1, '-A', adapter2, '--minimum-length', '25', '-j', '8', '-o', output1, '-p', output2, input1, input2]
	with open(statsout, 'w') as outfh:
		subprocess.call(command, stdout = outfh)