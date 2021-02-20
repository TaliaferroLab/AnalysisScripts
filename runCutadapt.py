#python3
import subprocess

samples = ['Ctrl.Rep1_1', 'Ctrl.Rep2_1', 'Ctrl.Rep3_1', 'KD.Rep1_1', 'KD.Rep2_1', 'KD.Rep3_1']


#adapter1 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter2 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'
adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
adapter3 = 'TGGTGGCTGGTGTGGCCAAGCTTCGATATCCGCATGCTA'
adapter4 = 'GGGAAAAAGATCTCAGTGGTATTTGTGAGCCAGCACTAGTGCGACCGCAAGAG'
#adapter5 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter6 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'


for idx, sample in enumerate(samples):
	print('Trimming {0}, sample {1} of {2}...'.format(sample, idx + 1, len(samples)))
	input1 = sample + '.fq.gz'
	input2 = input1.replace('_1.fq.gz', '_2.fq.gz')
	if 'Soma_FF' in input1:
		index = 13
	elif 'Neurite_FF' in input1:
		index = 16
	elif 'Soma_GFP' in input1:
		index = 14
	elif 'Neurite_GFP' in input1:
		index = 17
	
	#samplename = sample[:index]
	#samplename = samplename[:-2] + '_Rep' + samplename[-1]
	samplename = sample.split('_')[0]
	output1 = '{0}.R1.trimmed.fq.gz'.format(samplename)
	output2 = '{0}.R2.trimmed.fq.gz'.format(samplename)
	statsout = '{0}.cutadaptstats.txt'.format(samplename)
	if 'FF' in samplename:
		adapters = [adapter1, adapter2]
	elif 'GFP' in samplename:
		adapters = [adapter3, adapter4]
	#command = ['cutadapt', '-g', adapters[0], '-G', adapters[1], '--minimum-length', '75', '-j', '8',
	#'-o', output1, '-p', output2, input1, input2]
	command = ['cutadapt', '-a', adapter1, '-A', adapter2, '--minimum-length', '25', '-j', '8', '-o', output1, '-p', output2, input1, input2]
	with open(statsout, 'w') as outfh:
		subprocess.call(command, stdout = outfh)