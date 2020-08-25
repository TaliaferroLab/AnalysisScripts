#Must load modules: biopython, numpy, sqlalchemy, mysql-python, pycogent
#Holds genome sequence in memory.  Do not run on your laptop.

#This script is similar to motifconservation_v2.0.py, except that this one has
#removed the requirement for gffutils.  I was getting sqlite3 errors saying that the
#gff database was locked.  I had run into this previously when using gffutils on coyote,
#and had just instead run scripts on my laptop.  However, because this uses the ensembl
#server installed on coyote, it must be run there.

#Usage: python motifconservation.py <singlekmer / kmerlist> <kmer_of_interest or kmerenrichmentfile> <k> <simple/complex> <number_of_control_sequences> <gff_of_3UTRs> <genome.fa>

import sys
import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
#import gffutils
from cogent.db.ensembl import HostAccount, Genome, Compara
from sqlalchemy.exc import OperationalError as OE
from MySQLdb import OperationalError as mOE

#Use if kmer is complex. Returns shuffled versions of that kmer with matching CpG contents.
def shufflekmer(kmer, iterations):
    iterations = int(iterations)
    sys.stderr.write('Shuffling kmer {0} times and matching CpG count to create control kmers... \n'.format(iterations))
    controlsequences = []
    kmer = kmer.upper()
    CpGcount = kmer.count('CG')
    kmer = list(kmer)
    while len(controlsequences) < iterations:
        random.shuffle(kmer)
        shuffled = ''.join(kmer)
        if shuffled.count('CG') == CpGcount and shuffled != kmer:
            controlsequences.append(shuffled)

    sys.stderr.write('Done!\n')
    return controlsequences #kmerlist input for kmerhomology

#Use if kmer is not complex.  Returns random kmers with same GC and CpG contents.
def randomkmer(kmer, iterations):
    iterations = int(iterations)
    sys.stderr.write('Creating {0} random kmers with matched GC content and CpG counts for kmer {1}...\n'.format(iterations, kmer))
    controlsequences = []
    kmer = kmer.upper()
    CpGcount = kmer.count('CG')
    GCcontent = kmer.count('C') + kmer.count('G')
    kmersize = len(kmer)
    while len(controlsequences) < iterations:
        randomseq = ''.join(random.choice('ACGT') for x in range(kmersize))
        randomCpGcount = randomseq.count('CG')
        randomGCcontent = randomseq.count('C') + randomseq.count('G')
        if randomCpGcount == CpGcount and randomGCcontent == GCcontent and randomseq != kmer:
            controlsequences.append(randomseq)
            
    sys.stderr.write('Done!\n')
    return controlsequences #kmerlist input for kmerhomology

#Get sequences of UTRs of interest in fasta format.  Sequence name contains ID,chr,start,stop,strand
def getSequences(gff, sequence_file):
    sys.stderr.write('Indexing genome sequences...\n')
    seq_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, 'fasta'))
    sys.stderr.write('Indexed {0} sequences.\n'.format(len(seq_dict)))

    '''
    #Make gff database
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    seqs = {} #fasta format in dictionary {seqID:sequence}
    
    UTRs = db.features_of_type('oligo') #MAY NEED TO CHANGE!!!
   
    sys.stderr.write('Retrieving sequences from gff...\n')
    for UTR in UTRs:
        chrm = str(UTR.chrom)
        strand = str(UTR.strand)
        start = int(UTR.start)
        stop = int(UTR.stop)
        ID = str(UTR.id)
    '''

    seqs = {} #fasta format in dictionary {seqID:sequence}
    gfffh = open(gff, 'r')
    for line in gfffh:
        line = line.strip().split('\t')
        if line[1] == 'SE': #may need to change
            chrm = line[0]
            start = int(line[3])
            stop = int(line[4])
            strand = line[6]
            ID = line[8][3:].split(';')[0] #This is expecting this field to be ID=<...>
            
            
            if strand == '+':
                #MAY NEED TO REMOVE STOP CODON
                #seqID = ID + ';' + chrm + ';' + str(start) + ';' + str(stop-50) + ';' + strand #gonna take off last 50 nt
                seqID = ID + ';' + chrm + ';' + str(start) + ';' + str(stop) + ';' + strand #dont take off last 50 nt
            elif strand == '-':
                #MAY NEED TO REMOVE STOP CODON
                #seqID = ID + ';' + chrm + ';' + str(start+50) + ';' + str(stop) + ';' + strand #gonna take off last 50 nt
                seqID = ID + ';' + chrm + ';' + str(start) + ';' + str(stop) + ';' + strand #dont take off last 50 nt
            UTRsequence = ''

            if strand == '+':
                UTRsequence += seq_dict[chrm].seq[start-1:stop].upper() #MAY NEED TO REMOVE STOP CODONS
            elif strand == '-':
                UTRsequence += seq_dict[chrm].seq[start-1:stop].upper().reverse_complement() #MAY NEED TO REMOVE STOP CODONS
    
        #if len(UTRsequence) > 50: #Remove last 50 nt
            #seqs[seqID] = str(UTRsequence[:-50])

            seqs[seqID] = str(UTRsequence) #Dont remove last 50 nt

            #os.remove(db_fn)
    sys.stderr.write('Retrieved {0} UTR sequences over 50 nt.\n'.format(len(seqs)))
    gfffh.close()
    return seqs #fastdict input for kmerhomology

#Go kmer by kmer in each sequence and see if conserved.  Expects fasta file as a dictionary as input (output from getSequences).
def kmerhomology(k, kmerlist, fastadict):
    k = int(k)
    homologydict = {} #{kmer:[conserved_occurences, total occurences]}
    UTRcounter = 0
    account = HostAccount('sugarman', 'ensembl', 'ensembl')
    #account = None
    compara = Compara(['mouse', 'human'], Release=61, account=account)
    sqlalchemyfails = 0
    
    if k != len(kmerlist[0]):
        sys.stderr.write('Warning! Provided value of k does not match length of given kmer!')
    for UTR in fastadict:
        UTRcounter +=1
        if UTRcounter % 50 == 0:
            sys.stderr.write('Determining motif conservation in UTR {0} of {1}...\n'.format(UTRcounter, len(fastadict)))
        UTRsequence = fastadict[UTR]
        UTR = UTR.replace(';', '\t').split('\t')
        ID = UTR[0]
        chrm = UTR[1].replace('chr','') #change to ensembl style
        start = int(UTR[2])
        stop = int(UTR[3])
        strand = UTR[4]

        for i in range(len(UTRsequence) - k + 1):
            if strand == '+':
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = start + i
                mousekmerstop = start + i + k - 1
            elif strand == '-': #actually counting back from the end...and remember the last 50 and stop codons are removed!!
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = stop - i - k + 1
                mousekmerstop = stop - i
            if mousekmer in kmerlist:
                if homologydict.has_key(mousekmer) == False:
                    homologydict[mousekmer] = [0, 1] #if not in dictionary, initiate entry with 1 in total occurences
                elif homologydict.has_key(mousekmer):
                    homologydict[mousekmer][1] +=1 #if in dictionary, add one to total occurences
                    
                #stupid fucking sqlalchemy timeouts
                try:
                    for synt_region in compara.getSyntenicRegions(Species = 'mouse', CoordName = chrm, Start = mousekmerstart, End = mousekmerstop, Strand = strand, ensembl_coord = True, align_method = 'PECAN', align_clade = '19 amniota vertebrates Pecan'):
                        membs = synt_region.Members
                        if len(membs) == 2: #if there is no aligned human seq just skip it
                            completed = True
                            mouse = membs[0]
                            human = membs[1]
                            #strand sometimes seems to not get picked up by compara. If on minus strand, this seq may be rev comp of motif.
                            mouseseq = str(mouse.AlignedSeq) 
                            humanseq = str(human.AlignedSeq)
                            if mouseseq == humanseq:
                                homologydict[mousekmer][0] += 1 #add conserved occurence to dictionary
                        elif len(membs) != 2:
                            pass
                except (OE, mOE):
                    sys.stderr.write('Genome mysql error!!!')
                    sqlalchemyfails +=1
                    i = i-1 #try again
                    continue
            
    return homologydict, sqlalchemyfails
            
def summarize_homologydict(homologydict, kmer_of_interest):
    controlkmer_conservedfractions = []

    for kmer in homologydict:
        conservedfraction = (homologydict[kmer][0] / float(homologydict[kmer][1]))
        if kmer != kmer_of_interest:
            controlkmer_conservedfractions.append(conservedfraction)
        elif kmer == kmer_of_interest:
            kmer_conservedfraction = conservedfraction

    controlkmer_mean = np.mean(controlkmer_conservedfractions)
    controlkmer_std = np.std(controlkmer_conservedfractions)
    zscore = ((kmer_conservedfraction - controlkmer_mean) / float(controlkmer_std))
    sys.stderr.write('For kmer {0}, the fraction of kmer occurences that were conserved was {1} while for the control kmers it was {2} with a stdev of {3}.\n'.format(kmer_of_interest, kmer_conservedfraction, controlkmer_mean, controlkmer_std))
    sys.stderr.write('The z-score for kmer {0} is {1}.\n'.format(kmer_of_interest, zscore))

    return zscore

def summarize_kmerhomologydict(kmerhomologydict, kmer_of_interest, kmerlist):
    controlkmer_conservedfractions = []
    if kmerhomologydict.has_key(kmer_of_interest) == False: #if kmer not in UTRs
        kmer_conservedfraction = 'NA'
        controlkmer_mean = 'NA'
        controlkmer_std = 'NA'
        zscore = 'NA'
    elif kmerhomologydict.has_key(kmer_of_interest):
        for kmer in kmerhomologydict:
            if kmer != kmer_of_interest and kmer in kmerlist:
                conservedfraction = (kmerhomologydict[kmer][0] / float(kmerhomologydict[kmer][1]))
                controlkmer_conservedfractions.append(conservedfraction)
            elif kmer == kmer_of_interest and kmer in kmerlist:
                conservedfraction = (kmerhomologydict[kmer][0] / float(kmerhomologydict[kmer][1]))
                kmer_conservedfraction = conservedfraction

        controlkmer_mean = np.mean(controlkmer_conservedfractions)
        controlkmer_std = np.std(controlkmer_conservedfractions)
        if controlkmer_std != 0:
            zscore = ((kmer_conservedfraction - controlkmer_mean) / float(controlkmer_std))
        elif controlkmer_std == 0:
            zscore = 'NA_because_control_std_was_0'

    kmerstats = [str(kmer_conservedfraction), str(controlkmer_mean), str(controlkmer_std), str(zscore)]
        
    return kmerstats
        
def kmerhomologydict(k, fastadict): #for multiple kmers at once, homologydict of every kmer in every UTR sequence
    k = int(k)
    kmerhomologydict = {} #{kmer:[conserved_occurences, total_occurences]}
    UTRcounter = 0
    analyzedUTRs = 0
    account = HostAccount('sugarman', 'ensembl', 'ensembl')
    #account = None
    compara = Compara(['mouse', 'human'], Release=61, account=account)
    sqlalchemyfails = 0

    for UTR in fastadict:
        UTRcounter +=1
        UTRsequence = fastadict[UTR]
        UTR = UTR.replace(';', '\t').split('\t')
        ID = UTR[0]
        chrm = UTR[1].replace('chr','') #change to ensembl style
        start = int(UTR[2])
        stop = int(UTR[3])
        strand = UTR[4]
        if UTRcounter % 1 == 0:
            sys.stderr.write('Determining motif conservation in UTR {0}, number {1} of {2}, (interrogated {3} so far)...\n'.format(ID, UTRcounter, len(fastadict), analyzedUTRs))

        for i in range(len(UTRsequence) - k + 1):
            if strand == '+':
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = start + i
                mousekmerstop = start + i + k - 1
            elif strand == '-': #actually counting back from the end...and remember the last 50 and stop codons are removed!!
                mousekmer = UTRsequence[i:i+k]
                mousekmerstart = stop - i - k + 1
                mousekmerstop = stop - i
            if kmerhomologydict.has_key(mousekmer) == False:
                kmerhomologydict[mousekmer] = [0, 1] #if not in dictionary, initiate entry with 1 in total occurences
            elif kmerhomologydict.has_key(mousekmer):
                kmerhomologydict[mousekmer][1] +=1 #if in dictionary, add one to total occurences
                
            for synt_region in compara.getSyntenicRegions(Species = 'mouse', CoordName = chrm, Start = mousekmerstart, End = mousekmerstop, Strand = strand, ensembl_coord = True, align_method = 'PECAN', align_clade = '19 amniota vertebrates Pecan'):
                if mousekmerstart == start and strand == '+' or mousekmerstop == stop and strand == '-': #this will be true once per UTR
                    analyzedUTRs +=1
                membs = synt_region.Members
                if len(membs) == 2: #if there is no aligned human seq just skip it
                    completed = True
                    mouse = membs[0]
                    human = membs[1]
                    #strand sometimes seems to not get picked up by compara. If on minus strand, this seq may be rev comp of motif.
                    mouseseq = str(mouse.AlignedSeq) 
                    humanseq = str(human.AlignedSeq)
                    if mouseseq == humanseq:
                        kmerhomologydict[mousekmer][0] += 1 #add conserved occurence to dictionary
                elif len(membs) != 2:
                    pass
                
    sys.stderr.write('Analyzed {0} of {1} UTRs.\n'.format(analyzedUTRs, len(fastadict)))

    return kmerhomologydict, sqlalchemyfails

            
    

if __name__ == '__main__':
    if len(sys.argv) != 8:
        sys.stderr.write('You haven\'t entered the correct number arguments.  Exiting.\n')
        sys.exit()
    elif len(sys.argv) == 8:
        if sys.argv[1] == 'singlekmer':
            sys.stderr.write('Single kmer mode.\n')
            kmer_of_interest = sys.argv[2].upper().replace('U','T') #genome sequence has Ts not Us
            k = int(sys.argv[3])
            kmer_complexity = sys.argv[4]
            control_iterations = int(sys.argv[5])
            UTRgff = sys.argv[6]
            genomefasta = sys.argv[7]

            if kmer_complexity == 'simple':
                kmerlist = randomkmer(kmer_of_interest, control_iterations) + [kmer_of_interest]
            elif kmer_complexity == 'complex':
                kmerlist = shufflekmer(kmer_of_interest, control_iterations) + [kmer_of_interest]

            fastadict = getSequences(UTRgff, genomefasta)
            motifconservationdict = kmerhomology(k, kmerlist, fastadict)[0]
            summarize_homologydict(motifconservationdict, kmer_of_interest)

        elif sys.argv[1] == 'kmerlist': #expects output from kmerenrichment.py
            sqlalchemyfails = 0
            kmerlistfile = open(sys.argv[2], 'r')
            k = int(sys.argv[3])
            kmer_complexity = sys.argv[4]
            control_iterations = int(sys.argv[5])
            UTRgff = sys.argv[6]
            genomefasta = sys.argv[7]
            fastadict = getSequences(UTRgff, genomefasta)
            kmerhomologydict, fails = kmerhomologydict(k, fastadict)
            print kmerhomologydict
            sqlalchemyfails += fails
            fh = open('kmerconservations2.txt', 'w')
            fh.close()
            for line in kmerlistfile:
                line = line.strip().split('\t')
                if line[0] != 'kmer':
                    kmer_of_interest = line[0].upper().replace('U','T')
                    if kmer_of_interest == 'CGCGCG': #impossible to get control kmers for this
                        continue
                    kmerenrichment = str(line[3])
                    pvalue = str(line[4])
                    bhpvalue = str(line[5])
                    if kmer_complexity == 'simple':
                        kmerlist = randomkmer(kmer_of_interest, control_iterations) + [kmer_of_interest]
                    elif kmer_complexity == 'complex':
                        kmerlist = shufflekmer(kmer_of_interest, control_iterations) + [kmer_of_interest]

                    #Summarize with summarize_kmerhomologydict
                    kmerstats = summarize_kmerhomologydict(kmerhomologydict, kmer_of_interest, kmerlist)
                    outfh = open('kmerconservations2.txt', 'a')
                    outfh.write(kmer_of_interest.replace('T','U') + '\t' + line[1] + '\t' + line[2] + '\t' + kmerenrichment + '\t' + pvalue
+ '\t' + bhpvalue + '\t' + kmerstats[0] + '\t' + kmerstats[1] + '\t' + kmerstats[2] + '\t' + kmerstats[3] + '\n')
                    outfh.close()
            sys.stderr.write('There were {0} sqlalchemy timeouts.\n'.format(sqlalchemyfails))
            sys.stderr.write('All done!\n')
            kmerlistfile.close()
            
