import argparse
import subprocess
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from UTR3gfftosequence_seqdict import gfftosequence as gfftosequence

parser = argparse.ArgumentParser()
parser.add_argument('--events', type = str, help = 'List of MISO events to process.')
parser.add_argument('--name', type = str, help = 'Group name for events.')
parser.add_argument('--annotations', type = str, help = 'File containing MISO annotations for events.')
parser.add_argument('--stopcodongff', type = str, help = 'Gff file containing stop codons explicitly annotated as \'stop_codon\'. mm9.genes.chr.jess.gff works well.')
parser.add_argument('--genomefasta', type = str, help = 'Fasta file of genome.')
args = parser.parse_args()
events = args.events
name = args.name
annotations = args.annotations
stopcodons = args.stopcodongff
genomefasta = args.genomefasta


#Function that returns number of lines in a file
def getlinecount(filename):
    output = subprocess.check_output('wc ' + '-l ' + filename, shell='TRUE')
    linecount = output.strip().split(' ')[0]
    return linecount

#First, get MISO event annotations for the events using RetrieveMISOannotations.py
#Sorting of annotations (proximal/distal) is lost.
print 'Starting with {0} events.'.format(getlinecount(events))
print 'Retrieving annotations for events...'
os.system('python' + ' /Users/mtaliaferro/Scripts/RetrieveMISOAnnotations.py ' + events + ' ' + annotations + ' ' + name + '.gff3')
eventannotations = name + '.gff3'

#Now resort them to put the distal isoform as the first listed isoform
#This assumes all events are two-isoform events
print 'Sorting retrieved annotations...'
os.system('python' + ' /Users/mtaliaferro/Scripts/ReorderGff3v2.0.py ' + eventannotations + ' ALE ' + '> ' + name + '_sorted.gff3')
sortedannotations = name + '_sorted.gff3'

#Using awk, split annotations into distal (first isoform) and proximal (second isoform) events
print 'Splitting annotations into distal and proximal isoforms...'
os.system('grep' + ' mRNA ' + sortedannotations + ' | ' + 'awk' + ' \'NR%2==1\' ' + '> ' + name + '_distalisoforms.gff3')
os.system('grep' + ' mRNA ' + sortedannotations + ' | ' + 'awk' + ' \'NR%2==0\' ' + '> ' + name + '_proximalisoforms.gff3')
distalisoforms = name + '_distalisoforms.gff3'
proximalisoforms = name + '_proximalisoforms.gff3'
print '{0} distal isoforms.'.format(getlinecount(distalisoforms))
print '{0} proximal isoforms.'.format(getlinecount(proximalisoforms))

#Using awk, collapse events that are duplicate (contain non-unique chromosome AND start AND stop AND strand entries) to one instance.
print 'Collapsing duplicate events...'
os.system('awk \'{ if (a[$1,$4,$5,$7]++ ==0) print $0; }\'' + ' ' + distalisoforms + ' > ' + name + '_distalisoforms_collapsed.gff3')
os.system('awk \'{ if (a[$1,$4,$5,$7]++ ==0) print $0; }\'' + ' ' + proximalisoforms + ' > ' + name + '_proximalisoforms_collapsed.gff3')
distalisoforms = name + '_distalisoforms_collapsed.gff3'
proximalisoforms = name + '_proximalisoforms_collapsed.gff3'
print '{0} collapsed distal isoforms.'.format(getlinecount(distalisoforms))
print '{0} collapsed proximal isoforms.'.format(getlinecount(proximalisoforms))

#Using getUTRcoords.py and a reference gff (mm9.genes.chr.jess.gff works well) get regions of events that lie outside of stop codon.
print 'Using stop codon annotation to determine UTR portions of events...'
os.system('python' + ' /Users/mtaliaferro/Scripts/get3UTRcoords.py ' + distalisoforms + ' ' + stopcodons + ' ' + name + '_distalisoformsUTRs.gff3')
os.system('python' + ' /Users/mtaliaferro/Scripts/get3UTRcoords.py ' + proximalisoforms + ' ' + stopcodons + ' ' + name + '_proximalisoformsUTRs.gff3')
distalisoformUTRs = name + '_distalisoformsUTRs.gff3'
proximalisoformUTRs = name + '_proximalisoformsUTRs.gff3'
print 'Found UTRs for {0} of {1} distal isoforms.'.format(getlinecount(distalisoformUTRs), getlinecount(distalisoforms))
print 'Found UTRs for {0} of {1} proximal isoforms.'.format(getlinecount(proximalisoformUTRs), getlinecount(proximalisoforms))


#Retrieve sequences for UTRs using 3UTRgfftosequence_seqdict.py.  This also removes last 50 nt of UTR, so it requires that the UTR be at least 51 nt long.
print 'Retrieving sequences for UTRs and removing last 50 nt...'
seqdictionary = gfftosequence(distalisoformUTRs, genomefasta)
outfh = open(name + '_distalisoformsUTRs.fasta', 'w')
for ID in seqdictionary:
        outfh.write('>' + ID + '\n' + str(seqdictionary[ID]) + '\n')
outfh.close()
seqdictionary = gfftosequence(proximalisoformUTRs, genomefasta)
outfh = open(name + '_proximalisoformsUTRs.fasta', 'w')
for ID in seqdictionary:
        outfh.write('>' + ID + '\n' + str(seqdictionary[ID]) + '\n')
outfh.close()
distalisoformseqs = name + '_distalisoformsUTRs.fasta'
proximalisoformseqs = name + '_proximalisoformsUTRs.fasta'
distalseqcount = int(getlinecount(distalisoformseqs)) / 2
proximalseqcount = int(getlinecount(proximalisoformseqs)) / 2
print 'Found UTRs of at least 51 nt for {0} distal ALEs and {1} proximal ALEs.'.format(str(distalseqcount), str(proximalseqcount))

#Now require that for a UTR sequence to be reported, the proximal AND distal UTR sequences must have been found and passed the length filter of 50 nt.
os.system('python' + ' /Users/mtaliaferro/Scripts/requireBothExons.py ' + distalisoformseqs + ' ' + proximalisoformseqs + ' 0 ' + name + '_distalisoformsUTRs_filtered.fasta ' + name + '_proximalisoformsUTRs_filtered.fasta')
distalisoformseqs = name + '_distalisoformsUTRs_filtered.fasta'
proximalisoformseqs = name + '_proximalisoformsUTRs_filtered.fasta'
distalseqcount = int(getlinecount(distalisoformseqs)) / 2
proximalseqcount = int(getlinecount(proximalisoformseqs)) / 2
if distalseqcount != proximalseqcount:
    print 'ERROR: After filtering, number of reported distal and proximal UTR sequences should equal each other.'
    #sys.exit()
print 'Finally, reported {0} UTR sequences for ALE events in which both proximal and distal UTR sequences passed length filters.'.format(distalseqcount)

#Cleanup
if os.path.isfile(events.split('.txt')[0] + '.gff3.db-shm'):
    os.remove(events.split('.txt')[0] + '.gff3.db-shm')
if os.path.isfile(events.split('.txt')[0] + '.gff3.db-wal'):
    os.remove(events.split('.txt')[0] + '.gff3.db-wal')
