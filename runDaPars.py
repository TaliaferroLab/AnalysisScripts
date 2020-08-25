import gffutils
import os
import sys
import argparse
import subprocess

#This script takes a gff (it's setup to take one from Gencode) and returns
#the annotation in bed12 format for all transcripts in the gff as well as a file
#that relates all transcripts to their parent genes.  These files are necessary for analysis
#of alternative polyadenylation using DaPars.

#Convert all 'transcripts' in a gff to bed12 format.  This can be used as input for the first DaPars script.
def gff2bed12(gff, bedoutfile, nameoutfile):
	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	transcripts = db.features_of_type('transcript')

	transcriptcounter = 0
	with open(bedoutfile, 'w') as bedoutfh, open(nameoutfile, 'w') as nameoutfh:
		for transcript in transcripts:
			transcriptcounter +=1
			if transcriptcounter % 1000 == 0:
				print 'Converting transcript {0}...'.format(transcriptcounter)
			bed12line = gffutils.FeatureDB.bed12(db, transcript)
			bed12line = bed12line.split('\t')
			bed12line[3] = bed12line[3].split('.')[0]
			bedoutfh.write(('\t').join(bed12line) + '\n')
			txid = transcript.id.split('.')[0]
			for gene in db.parents(transcript, featuretype = 'gene'):
				geneid = gene.id.split('.')[0]
			nameoutfh.write(txid + '\t' + geneid + '\n')

def bam2bedgraph(bam):
	outfilename = ('.').join(os.path.basename(bam).split('.')[:-1]) + '.bedgraph'
	command = ['bedtools', 'genomecov', '-bg', '-split', '-ibam', bam]
	with open(outfilename, 'w') as outfh:
		subprocess.call(command, stdout=outfh)

def makeconfigfile(UTR3bed, cond1bams, cond2bams, expname):
	#UTR3bed = the output from the first DaPars script
	#cond1bams = comma separated list of bams for condition1
	#cond2bams = comma separated list of bams for condition2
	#expname = Name for this experiment

	#Have to take bam names and change them to the names of the bedgraphs you get out after running bam2bedgraph
	cond1bams = [('.').join(os.path.basename(bam).split('.')[:-1]) + '.bedgraph' for bam in cond1bams]
	cond2bams = [('.').join(os.path.basename(bam).split('.')[:-1]) + '.bedgraph' for bam in cond2bams]

	configfn = 'DaParsConfig_' + expname + '.txt'

	with open(configfn, 'w') as outfh:
		outfh.write('Annotated_3UTR={0}'.format(UTR3bed) + '\n')
		outfh.write('Group1_Tophat_aligned_Wig={0}'.format((',').join(cond1bams)) + '\n')
		outfh.write('Group2_Tophat_aligned_Wig={0}'.format((',').join(cond2bams)) + '\n')
		outfh.write('Output_directory=DaParsOutput/' + '\n')
		outfh.write('Output_result_file={0}'.format(expname) + '\n')
		outfh.write('Num_least_in_group1=1' + '\n')
		outfh.write('Num_least_in_group2=1' + '\n')
		outfh.write('Coverage_cutoff=30' + '\n')
		outfh.write('FDR_cutoff=0.05' + '\n')
		outfh.write('PDUI_cutoff=0.10' + '\n')
		outfh.write('Fold_change_cutoff=0' + '\n')

def runDaPars(UTR3bed, cond1bams, cond2bams, expname):
	#for each bam, run bam2bedgraph
	#make config file
	#run dapars script
	for bam in cond1bams:
		print 'Converting {0} to bedgraph...'.format(os.path.basename(bam))
		bam2bedgraph(bam)
	for bam in cond2bams:
		print 'Converting {0} to bedgraph...'.format(os.path.basename(bam))
		bam2bedgraph(bam)

	makeconfigfile(UTR3bed, cond1bams, cond2bams, expname)
	configfn = 'DaParsConfig_' + expname + '.txt'

	command = ['python', '/vol3/home/taliaferro/Tools/dapars/src/DaPars_main.py', configfn]
	print 'Running DaPars...'
	subprocess.call(command)
	print 'All done!'


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', type = str, help = 'Mode.', choices = ['makeUTRs', 'runDaPars'])
	parser.add_argument('--gff', type = str, help = 'Annotation in gff format.')
	parser.add_argument('--bed12', type = str, help = 'Bed12 output.')
	parser.add_argument('--tx2gene', type = str, help = 'Transcript/gene relations output.')
	parser.add_argument('--cond1bams', type = str, help = 'Comma separated list of sorted bam files from Condition 1.')
	parser.add_argument('--cond2bams', type = str, help = 'Comma spearated list of sorted bam files from Condition 2.')
	parser.add_argument('--UTR3bed', type = str, help = 'Output from first DaPars script.')
	parser.add_argument('--expname', type = str, help = 'Name of this experiment.')
	args = parser.parse_args()

	if args.mode == 'makeUTRs':
		gff2bed12(args.gff, args.bed12, args.tx2gene)
	elif args.mode == 'runDaPars':
		runDaPars(args.UTR3bed, args.cond1bams.split(','), args.cond2bams.split(','), args.expname)
