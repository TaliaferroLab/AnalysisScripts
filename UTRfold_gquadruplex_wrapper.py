import subprocess
import os

for fasta in os.listdir('.'):
	if os.path.basename(fasta).endswith('fa'):
		print fasta
		name = os.path.basename(fasta).split('.')[0] + '.' + os.path.basename(fasta).split('.')[1]
		output = name + '.gquadout.txt'
		windowsize = '80'
		slidesize = '10'
		print 'Folding {0}...'.format(name)
		command = ['python', '/home/jmtali/Scripts/UTRfold_gquadruplex.py', '--fasta', fasta, '--output', output, '--windowsize', windowsize, '--slidesize', slidesize]
		subprocess.call(command)
