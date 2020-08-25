#Calculates average phastcons scores for a window.  Then slides that window and recalculates.  Takes the max score for each UTR as the score for that UTR.  Phastcons scores are expected to be in a bed and tabix indexed.  UTRs are expected to be in a gff.

#Usage: python phastconswindow.py ALEDistal0.1AxonDISTALUTRs_filtered.gff3 sorted.phastcons.mm9.bed.gz outfile.txt



import sys
import os
import gffutils
import tabix
import pysam

def phastconswindow(gff, phastconsbed):
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    phastconstabix = pysam.Tabixfile(phastconsbed)

    maxwindowscores = {} # {UTR_ID : max score in any window}
    allwindowscores = {} # {UTR_ID : [all window scores]}
    interrogated_3UTRs = 0
    UTRs_with_scores = 0
    three_prime_UTRs = db.features_of_type('3\'UTR')

    for UTR in three_prime_UTRs:
        interrogated_3UTRs +=1
        if interrogated_3UTRs % 50 == 0:
            sys.stderr.write('Interrogating UTR {0}... \n'.format(interrogated_3UTRs))
        currentlocation = int(UTR.start)
        averagewindowscores = []
        UTRlength = int(UTR.stop - UTR.start)
        if UTRlength >= 200:
            while currentlocation + 200 <= UTR.stop:
                windowscores = [] #score of every nt in the window
                window = [str(UTR.chrom), currentlocation, currentlocation + 50] #50 nt window size
                for bed in phastconstabix.fetch(window[0], window[1], window[2], parser = pysam.asBed()):
                    windowscores.append(float(bed.name))
                if len(windowscores) > 0: #sometimes even though there should be a window, there are no phastcons scores in that window
                    averagewindowscores.append(sum(windowscores) / float(len(windowscores)))
                
                currentlocation +=1 #1 nt step size

        if len(averagewindowscores) > 0:
            UTRs_with_scores +=1
            maxwindowscores[UTR.id] = max(averagewindowscores)
            allwindowscores[str(UTR.id)] = averagewindowscores

    sys.stderr.write('Searched {0} UTRs and got PhastCons scores for {1} of them. \n'.format(interrogated_3UTRs, UTRs_with_scores))

    os.remove(db_fn)
    #return maxwindowscores
    return allwindowscores

if __name__ == '__main__':
    outfh = open(sys.argv[3], 'w')
    outfh.write('Gene' + '\t' + 'Phastcons_windowscore' + '\n')
    maxwindowscores = phastconswindow(sys.argv[1], sys.argv[2])
    for UTR in maxwindowscores:
        outfh.write(str(UTR) + '\t' + str(maxwindowscores[UTR]) + '\n')
    outfh.close()
                
