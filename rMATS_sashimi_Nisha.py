import os
import subprocess
import sys
import pandas as pd
import argparse
import shutil

def makebamentries(bamdirectory, sample1, sample2):
	sample1bams = []
	sample2bams = []
	for rep in ['03', 'C1', 'C3']:
		sample1bam = os.path.join(bamdirectory, '{0}Aligned.sortedByCoord.out.bam'.format(rep))
		sample1bams.append(sample1bam)
	for rep in ['FXS', 'FXSB2', 'SO10C3']:
		sample2bam = os.path.join(bamdirectory, '{0}Aligned.sortedByCoord.out.bam'.format(rep))
		sample2bams.append(sample2bam)

	sample1bams = (',').join(sample1bams)
	sample2bams = (',').join(sample2bams)

	print sample1bams, sample2bams
	return sample1bams, sample2bams

def makeoutdir(sample1, sample2):
	if os.path.exists(os.path.join('.', 'SashimiPlots')) == False:
		os.mkdir(os.path.join('.', 'SashimiPlots'))
	os.mkdir(os.path.join('.', 'SashimiPlots', '{0}vs{1}'.format(sample1, sample2)))
	outdir = os.path.join('.', 'SashimiPlots', '{0}vs{1}'.format(sample1, sample2))
	return outdir

def makeeventsfile(events, rMATSdir, sample1, sample2):
	rmatsfile = os.path.join(rMATSdir, 'SE.MATS.JC.txt')
	df = pd.read_table(rmatsfile)
	eventdfs = []
	for event in events.split(','):
		event = int(event)
		eventdf = df.loc[(df['ID'] == event)]
		eventdfs.append(eventdf)
	outdf = pd.concat(eventdfs)
	with open('sashimievents.txt', 'w') as outfh:
		outdf.to_csv(outfh, sep = '\t', index = False)

def runsashimi(bamdirectory, events, rMATSdir, sample1, sample2, colors):
	sample1bams, sample2bams = makebamentries(bamdirectory, sample1, sample2)
	makeeventsfile(events, rMATSdir, sample1, sample2)
	outdir = makeoutdir(sample1, sample2)
	colors = colors.split(',')
	colors = ['#' + color for color in colors]
	colors = (',').join(colors)

	command = ['rmats2sashimiplot', '--b1', sample1bams, '--b2', sample2bams, '-t', 'SE', '-e', 'sashimievents.txt', '--l1', sample1, '--l2', sample2, '-o', outdir, '--color', colors]
	subprocess.call(command)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--bamdirectory', type = str, help = 'Directory containing bam files.')
	parser.add_argument('--events', type = str, help = 'Event ids separated by commas.  E.g. chrX:-:13276434-13276547:13282165-13282229:13285202-13285632')
	parser.add_argument('--rMATSdir', type = str, help = 'Directory that contains all rMATS comparison subdirectories.')
	parser.add_argument('--sample1', type = str, help = 'Sample 1 name.')
	parser.add_argument('--sample2', type = str, help = 'Sample 2 name.')
	parser.add_argument('--colors', type = str, help = 'Hex colors separated by commas without leading #.')
	args = parser.parse_args()

	runsashimi(args.bamdirectory, args.events, args.rMATSdir, args.sample1, args.sample2, args.colors)
	for subdir in os.listdir('.'):
		if os.path.isdir(os.path.join('.', subdir)) and subdir.startswith('Sashimi_index'):
			shutil.rmtree(subdir)
