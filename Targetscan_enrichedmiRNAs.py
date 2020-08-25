#This script is hijacked from targetscan_parsecontextscores.py.  It asks which miRNAs
#are enriched for having sites in a particular sequence set.  Actually, more precisely,
#it just gives the density of sites for each miRNA.  Number of sites for a miRNA / total sequence 
#search space.

import os
import gffutils
import argparse
from numpy import mean as mean
from numpy import median as median

def parsecontextscores(csfile, gff, featurename):
    #Make dictionary of this form:
    # {UTRname : [[UTRlength], [names of all miRNAs that have sites in that UTR]]}
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

    #Now get miRNA names
    csfilehandle = open(csfile, 'r')
    for line in csfilehandle:
        line = line.strip().split('\t')
        if line[0] != 'Gene ID': #skip header line
            featureid = line[0].split(';')[0] #Remove Parent=...
            species = line[1]
            miRNAname = line[2]
            if species == '10090': #this is mouse; for other species, change this number
                if featureid not in CSdict:
                    CSdict[featureid] = [[lengthdict[featureid]], [miRNAname]]
                elif featureid in CSdict:
                    CSdict[featureid][1].append(miRNAname)

    csfilehandle.close()
            
    return CSdict

def parseCSdict(CSdict):
    #CSdict = {UTRname : [[UTRlength], [names of all miRNAs that have sites in that UTR]]}
    miRNAsites = {} #{miRNA : number of sites}
    miRNAdensities = {} #{miRNA : density of sites}

    totalsequencelength = 0

    for UTR in CSdict:
        totalsequencelength += int(CSdict[UTR][0][0])

    print 'The total sequence search space was {0} nt'.format(totalsequencelength)

    #Count miRNA occurences
    for UTR in CSdict:
        miRNAs = CSdict[UTR][1]
        for miRNA in miRNAs:
            if miRNA not in miRNAsites:
                miRNAsites[miRNA] = 1
            elif miRNA in miRNAsites:
                miRNAsites[miRNA] +=1

    for miRNA in miRNAsites:
        miRNAdensities[miRNA] = miRNAsites[miRNA] / float(totalsequencelength)

    return miRNAdensities
                
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csfile', type = str, help = 'Targetscan_60_context_scores.pl output.')
    parser.add_argument('--gff', type = str, help = 'Gff of regions that targetscan looked through.')
    parser.add_argument('--featurename', type = str, help = 'Feature category in gff file (3rd field of gff)')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    CSdict = parsecontextscores(args.csfile, args.gff, args.featurename)
    miRNAdensities = parseCSdict(CSdict)
    outfh = open(args.outfile, 'w')
    outfh.write('miRNA' + '\t' + 'density' + '\n')
    for entry in miRNAdensities:
        outfh.write(entry + '\t' + str(miRNAdensities[entry]) + '\n')

    outfh.close()
    
    
