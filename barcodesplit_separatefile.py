from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import izip
import gzip
import sys

def splitpairedfastq_separatebarcodefile(forwardfastq, reversefastq, barcodefastq):
	forfh = gzip.open(forwardfastq)
	revfh = gzip.open(reversefastq)
	barcodefh = gzip.open(barcodefastq)
	readcounter = 0

	for ((fortitle, forseq, forqual), (revtitle, revseq, revqual), (barcodetitle, barcodeseq, barcodequal)) in izip(FastqGeneralIterator(forfh), FastqGeneralIterator(revfh), FastqGeneralIterator(barcodefh)):
		barcode = barcodeseq
		readcounter +=1
		if readcounter % 1000000 == 0:
			print 'Analyzing read {0}...'.format(readcounter)

		if sum(ii == jj for ii, jj in izip('GACGAC', barcode)) >=5:
			with open('Input1_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Input1_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('TAATCG', barcode)) >=5:
			with open('NoProt1_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('NoProt1_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('TACAGC', barcode)) >=5:
			with open('Pulldown1_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Pulldown1_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('CACGAT', barcode)) >=5:
			with open('Input2_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Input2_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('CACTCA', barcode)) >=5:
			with open('NoProt2_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('NoProt2_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('CAGGCG', barcode)) >=5:
			with open('Pulldown2_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Pulldown2_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		if sum(ii == jj for ii, jj in izip('CATGGC', barcode)) >=6:
			with open('Input3_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Input3_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('CATTTT', barcode)) >=5:
			with open('NoProt3_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('NoProt3_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

		elif sum(ii == jj for ii, jj in izip('CCAACA', barcode)) >=5:
			with open('Pulldown3_1.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (fortitle, forseq, forqual))
			with open('Pulldown3_2.fastq', 'a') as outfile:
				outfile.write('@%s\n%s\n+\n%s\n' % (revtitle, revseq, revqual))

	forfh.close()
	revfh.close()
	barcodefh.close()

splitpairedfastq_separatebarcodefile(sys.argv[1], sys.argv[2], sys.argv[3])