#python3
import subprocess

samples = []
celllines = ['CAD', 'N2A']
comps = ['Neurite', 'Soma']
dox = ['MinusDox', 'PlusDox']
reps = ['Rep1', 'Rep2', 'Rep3']

for c in celllines:
	for co in comps:
		for d in dox:
			for r in reps:
				samp = c + '.' + co + '.' + d + '.' + r
				samples.append(samp)


#adapter1 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter2 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'
adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGA'
#adapter3 = 'AAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGCGCTTTTTTTTTT'
#adapter4 = 'AAAAAAAAAAGCGCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTT'
#adapter5 = 'GGCGGAAAGATCGCCGTGTAAGTTTGCTTCGATATCCGCATGCTA'
#adapter6 = 'CCACTTATTTCATGGATACTTGGAATGGTTTCACTAGTGCGACCGCAAGAG'


for idx, sample in enumerate(samples):
	print('Trimming {0}, sample {1} of {2}...'.format(sample, idx + 1, len(samples)))
	input1 = '{0}_1.fq.gz'.format(sample)
	input2 = '{0}_2.fq.gz'.format(sample)
	output1 = '{0}.R1.trimmed.fq.gz'.format(sample[:-4])
	output2 = '{0}.R2.trimmed.fq.gz'.format(sample[:-4])
	statsout = '{0}.cutadaptstats.txt'.format(sample)
	#command = ['cutadapt', '-g', adapter1, '-G', adapter2, '--minimum-length', '30', '-j', '8',
	#'-o', output1, '-p', output2, input1, input2]
	command = ['cutadapt', '-a', adapter1, '-A', adapter2, '--minimum-length', '25', '-j', '8', '-o', output1, '-p', output2, input1, input2]
	with open(statsout, 'w') as outfh:
		subprocess.call(command, stdout = outfh)