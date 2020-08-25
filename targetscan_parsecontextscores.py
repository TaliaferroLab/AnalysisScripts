#Takes the output of targetscan_60_context_scores.pl and a gff of the regions of interest
#and returns a file containing the length of each region, the average all of the context scores
#for every miRNA that targets that region, the sum of all of the context scores for every miRNA
#that targets that region, and that sum normalized by the length of the region.

import os
import gffutils
import argparse
from numpy import mean as mean
from numpy import median as median

def parsecontextscores(csfile, gff, featurename):
    #Make dictionary of this form:
    # {UTRname : [[UTRlength], [context scores for all miRNA binding sites]]}
    #csfile = output of targetscan_60_context_scores.pl
    #gff = gff file of regions of interest
    #featurename = feature category in gff file (3rd field)
    lengthdict = {}
    CSdict = {}

    #First need to get lengths
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False: #if database doesn't exist, create it
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    features = db.features_of_type(featurename)

    for feature in features:
        featureid = feature.id
        featurelength = feature.stop - feature.start
        lengthdict[featureid] = featurelength

    os.remove(db_fn)

    #Now get context scores
    csfilehandle = open(csfile, 'r')
    for line in csfilehandle:
        line = line.strip().split('\t')
        if line[0] != 'Gene ID': #skip header line
            featureid = line[0].split(';')[0] #Remove Parent=...
            species = line[1]
            contextscore = line[11]
            if species == '10090' and contextscore != 'too_close': #this is mouse; for other species, change this number
                if featureid not in CSdict:
                    CSdict[featureid] = [[lengthdict[featureid]], [float(contextscore)]]
                elif featureid in CSdict:
                    CSdict[featureid][1].append(float(contextscore))

    csfilehandle.close()
            
    return CSdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csfile', type = str, help = 'Targetscan_60_context_scores.pl output.')
    parser.add_argument('--gff', type = str, help = 'Gff of regions that targetscan looked through.')
    parser.add_argument('--featurename', type = str, help = 'Feature category in gff file (3rd field of gff)')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    CSdict = parsecontextscores(args.csfile, args.gff, args.featurename)
    outfh = open(args.outfile, 'w')
    outfh.write('Sequence' + '\t' + 'Length' + '\t' + 'Numberofsites' + '\t' + 'Mean_context_score' + '\t' + 'Median_context_score' + '\t' + 'Total_context_score' + '\t' 'Normalized_score' + '\n')
    for entry in CSdict:
        length = CSdict[entry][0][0]
        MeanContextScore = float(sum(CSdict[entry][1]) / len(CSdict[entry][1]))
        MedianContextScore = median(CSdict[entry][1])
        TotalContextScore = sum(CSdict[entry][1])
        Numberofsites = len(CSdict[entry][1])
        NormalizedScore = float(TotalContextScore / length)
        outfh.write(entry + '\t' + str(length) + '\t' + str(Numberofsites) + '\t' + str(MeanContextScore) + '\t' + str(MedianContextScore) + '\t' + str(TotalContextScore) + '\t' + str(NormalizedScore) + '\n')

    outfh.close()
    
    
