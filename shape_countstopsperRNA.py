import argparse
from numpy import mean
from numpy import median

#Get the number of treated and untreated stops in each RNA. These are taken from reactivities.out

def stopsperRNA(reactivities):
    reactivitiesfh = open(reactivities, 'r')
    stopdict = {} # {seqname : [treatedstops, untreatedstops]
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] == 'sequence': #skip header
            continue
        seqname = line[0]
        treatedstops = int(line[4])
        untreatedstops = int(line[5])
        if seqname not in stopdict:
            stopdict[seqname] = [treatedstops, untreatedstops]
        elif seqname in stopdict:
            stopdict[seqname][0] += treatedstops
            stopdict[seqname][1] += untreatedstops

    reactivitiesfh.close()

    return stopdict

def RTstoppos(reactivities):
    reactivitiesfh = open(reactivities, 'r')
    stopposdict = {} # {stoppos : [treatedstops, untreatedstops]
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] == 'sequence': #skip header
            continue
        stoppos = line[2]
        treatedstops = int(line[4])
        untreatedstops = int(line[5])
        if stoppos not in stopposdict:
            stopposdict[stoppos] = [treatedstops, untreatedstops]
        elif stoppos in stopposdict:
            stopposdict[stoppos][0] += treatedstops
            stopposdict[stoppos][1] += untreatedstops

    reactivitiesfh.close()

    return stopposdict

def meanreactperpos(normalizedreactivities):
    reactivitiesfh = open(normalizedreactivities, 'r')
    reactdict = {} # {pos : [list of reactivities]}
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        if line[0] == 'sequence': #skip header
            continue
        stoppos = line[2]
        if int(stoppos) > 0:
            normalized_theta = float(line[8])
            if stoppos in reactdict:
                reactdict[stoppos].append(normalized_theta)
            elif stoppos not in reactdict:
                reactdict[stoppos] = [normalized_theta]

    reactivitiesfh.close()

    return reactdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reactivities', type = str, required = True, help = 'reactivities.out from spats run.')
    parser.add_argument('--treatedreads', type = float, required = False, help = 'Number of treated reads. Used for normalizing.')
    parser.add_argument('--untreatedreads', type = float, required = False, help = 'Number of untreated reads. Used for normalizing.')
    parser.add_argument('--mode', type = str, required = True, choices =['stopsperRNA', 'RTstoppos', 'meanthetaperpos'], help = 'Do you want the number of stops per RNA, where those stops are, or the mean reactivity at each position?')
    parser.add_argument('--outfile', type = str, required = True, help = 'Output file.')
    args = parser.parse_args()

    if args.mode == 'stopsperRNA':
        stopdict = stopsperRNA(args.reactivities)
        outfh = open(args.outfile, 'w')
        outfh.write(('\t').join(['RNA','treatedstops','untreatedstops','normalizedtreated','normalizeduntreated']) + '\n')
        for RNA in stopdict:
            outfh.write(RNA + '\t' + str(stopdict[RNA][0]) + '\t' + str(stopdict[RNA][1]) + '\t' + str((stopdict[RNA][0] / args.treatedreads) * 1000000) + '\t' + 
                        str((stopdict[RNA][1] / args.untreatedreads) * 1000000) + '\n')

        outfh.close()

    elif args.mode == 'RTstoppos':
        stopposdict = RTstoppos(args.reactivities)
        outfh = open(args.outfile, 'w')
        outfh.write(('\t').join(['RTstoppos','treatedstops','untreatedstops','normalizedtreated','normalizeduntreated']) + '\n')
        for stoppos in sorted(stopposdict):
            outfh.write(stoppos + '\t' + str(stopposdict[stoppos][0]) + '\t' + str(stopposdict[stoppos][1]) + '\t' + str((stopposdict[stoppos][0] / args.treatedreads) * 1000000) + '\t' + 
                        str((stopposdict[stoppos][1] / args.untreatedreads) * 1000000) + '\n')

        outfh.close()

    elif args.mode == 'meanthetaperpos':
        print 'Make sure this is a normalized reactivity file!'
        reactdict = meanreactperpos(args.reactivities)
        outfh = open(args.outfile, 'w')
        outfh.write('RTstoppos' + '\t' + 'mean_normalized_theta' + '\n')
        for nt in sorted(reactdict, key = int):
            outfh.write(nt + '\t' + str(mean(reactdict[nt])) + '\n')
        outfh.close()

        
