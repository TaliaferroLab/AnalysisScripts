import re
from Bio import SeqIO
import sys

def getnumberofGquads(seq):
	seqlen = len(seq)
	#Look for all possible G-quadruplex forming seqs WGGA-(N0-6)-WGGA-(N0-6)-WGGA-(N0-6)-WGGA
	#Allow overlapping matches to count individually
	motif1matches = re.findall(r'(?=([AU]GGA(.{0,6})[AU]GGA(.{0,6})[AU]GGA(.{0,6})[AU]GGA))', seq)
	#motif1matches = re.findall(r'(?=([AU]GG(.{0,7})[AU]GG(.{0,7})[AU]GG(.{0,7})[AU]GG))', seq)
	return len(motif1matches)


with open(sys.argv[2], 'w') as f:
	f.write(('\t').join(['gene', 'WGGAgquadcount', 'seqlen', 'WGGAgquaddens']) + '\n')

for record in SeqIO.parse(sys.argv[1], 'fasta'):
	ggacount = str(record.seq).count('GGA')
	gquadcount = getnumberofGquads(str(record.seq.transcribe()))
	seqlen = len(str(record.seq))
	if seqlen == 0:
		continue
	if ggacount > 0:
		gquadfrac = gquadcount / float(ggacount)
	elif ggacount == 0:
		gquadfrac = 'NA'

	#with open(sys.argv[2], 'a') as f:
		#f.write(('\t').join([str(gquadcount), str(ggacount), str(seqlen), str(gquadcount / float(seqlen)), str(gquadfrac), sys.argv[3], sys.argv[4]]) + '\n')

	with open(sys.argv[2], 'a') as f:
		f.write(('\t').join([str(record.id.split('.')[0]), str(gquadcount), str(seqlen), str(gquadcount / float(seqlen))]) + '\n')
