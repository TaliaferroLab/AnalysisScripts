#python3
import os
import subprocess


samples = ['CAD_Neurite_Rep1_S68', 'CAD_Neurite_Rep2_S69', 'CAD_Neurite_Rep3_S70', 'CAD_Neurite_Rep4_S71', 'CAD_Neurite_Rep5_S72', 
'CAD_Soma_Rep1_S63', 'CAD_Soma_Rep2_S64', 'CAD_Soma_Rep3_S65', 'CAD_Soma_Rep4_S66', 'CAD_Soma_Rep5_S67', 
'N2A_Neurite_Rep1_S58', 'N2A_Neurite_Rep2_S59', 'N2A_Neurite_Rep3_S60', 'N2A_Neurite_Rep4_S61', 'N2A_Neurite_Rep5_S62', 
'N2A_Soma_Rep1_S53', 'N2A_Soma_Rep2_S54', 'N2A_Soma_Rep3_S55', 'N2A_Soma_Rep4_S56', 'N2A_Soma_Rep5_S57']

samplenames = ['CADNeuriteRep1', 'CADNeuriteRep2', 'CADNeuriteRep3', 'CADNeuriteRep4', 'CADNeuriteRep5',
'CADSomaRep1', 'CADSomaRep2', 'CADSomaRep3', 'CADSomaRep4', 'CADSomaRep5',
'N2ANeuriteRep1', 'N2ANeuriteRep2', 'N2ANeuriteRep3', 'N2ANeuriteRep4', 'N2ANeuriteRep5',
'N2ASomaRep1', 'N2ASomaRep2', 'N2ASomaRep3', 'N2ASomaRep4', 'N2ASomaRep5']

readdir = '/beevol/home/taliaferro/data/cisElementScreen/Fractionation/EqualRNAamount/RawReads/trimmed'
indexfile = '/beevol/home/taliaferro/data/cisElementScreen/Fractionation/EqualRNAamount/Alignments/Bowtie2Index/mm10oligos'
outputdir = '/beevol/home/taliaferro/data/cisElementScreen/Fractionation/EqualRNAamount/Alignments'



for idx, sample in enumerate(samples):
	print('Aligning {0}, sample {1} of {2}...'.format(sample, idx + 1, len(samples)))
	forreads = os.path.join(readdir, '{0}.R1.trimmed.fq.gz'.format(sample))
	revreads = os.path.join(readdir, '{0}.R2.trimmed.fq.gz'.format(sample))
	samname = os.path.join(outputdir, samplenames[idx] + '.sam')
	statsout = os.path.join(outputdir, '{0}.bowtiestats.txt'.format(samplenames[idx]))
	command = ['bowtie2', '-q', '--end-to-end', '--fr', '--no-discordant', '--no-unal', '-p', '8', '-x', indexfile, '-1', forreads, '-2', revreads, '-S', samname]
	with open(statsout, 'w') as outfh:
		subprocess.call(command, stderr = outfh)