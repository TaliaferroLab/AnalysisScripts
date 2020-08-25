#Given a fasta of CDS sequences, output a fasta of AA sequences.

from Bio import SeqIO
import sys

def translatefasta(fasta, output):
	outfh = open(output, 'w')

	for CDS in SeqIO.parse(fasta, 'fasta'):
		AAseq = str(CDS.seq.translate())
		if len(CDS.seq) % 3 != 0:
			print 'WARNING: The length of the CDS {0} is not divisible by 3!! Skipping!!'.format(CDS.id)
			if AAseq[-1] != '*':
				print 'WARNING: CDS {0} does not end with a stop codon!! Skipping!!'.format(CDS.id)
			continue

		if AAseq[-1] != '*':
			print 'WARNING: CDS {0} does not end with a stop codon!! Skipping!!'.format(CDS.id)
			continue

		#Remove stop codon mark (*)
		AAseq = AAseq[:-1]

		outfh.write('>' + CDS.id + '\n' + AAseq + '\n')

	outfh.close()


translatefasta(sys.argv[1], sys.argv[2])