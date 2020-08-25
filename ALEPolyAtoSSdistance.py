import gffutils
import argparse
import os
import numpy

def getdistances(ALEevents, ALEannotation):
    distances = {} #{event : distance}
    distancelist = [] 
    eventsofinterest = []
    foundevents = 0
    eventsofinterestfh = open(ALEevents, 'r')
    for line in eventsofinterestfh:
        line = line.strip()
        eventsofinterest.append(line)
    eventsofinterestfh.close()
    proximalisoforms = [] 
    distalisoforms = [] 

    gff_fn = ALEannotation
    db_fn = os.path.basename(gff_fn) + '.db'

    print 'Indexing annotation...'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, force = False, verbose = False)
    db = gffutils.FeatureDB(db_fn)
    events = db.features_of_type('gene')

    for event in events:
        if event.id in eventsofinterest:
            foundevents +=1 #this event exists in the annotation
            if event.strand == '+':
                isoformcounter = 1 
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'end'): #order by end of isoform, this puts proximal isoform first
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter == 1:
                        proximalisoform = [isoform.chrom, 'ALEprox3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid]
                        isoformcounter +=1
                    elif isoformcounter == 2:
                        distalisoform = [isoform.chrom, 'ALEdist3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid]
                proxisoformpolyA = int(proximalisoform[4])
                distisoform3ss = int(distalisoform[3])

            elif event.strand == '-':
                isoformcounter = 1 
                #order by end of isoform, this puts proximal isoform first
                for isoform in db.children(event, featuretype = 'mRNA', order_by = 'start', reverse = True): 
                    for parent in db.parents(isoform, featuretype = 'gene'):
                        parentid = parent.id
                    if isoformcounter == 1:
                        proximalisoform = [isoform.chrom, 'ALEprox3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid]
                        isoformcounter +=1
                    elif isoformcounter == 2:
                        distalisoform = [isoform.chrom, 'ALEdist3ss', 'mRNA', str(isoform.start), str(isoform.end), '.', isoform.strand,
                                                 '.', isoform.id + ';Parent=' + parentid]

                proxisoformpolyA = int(proximalisoform[3])
                distisoform3ss = int(distalisoform[4])
                
            distance = abs(proxisoformpolyA - distisoform3ss)
            distancelist.append(distance)
            distances[parentid] = distance
    
    print 'Found {0} of {1} provided events in the annotation.'.format(foundevents, len(eventsofinterest))
    
    print 'The mean polyA to 3\' ss distance for these events is {0} +\- {1}'.format(numpy.mean(distancelist), numpy.std(distancelist))
    

    #Cleanup
    os.remove(db_fn)

    return distances

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type = str, help = 'List of ALE events to consider.', required = True)
    parser.add_argument('--annotations', type = str, help = 'ALE annotations containing events in MISO gff format.', required = True)
    parser.add_argument('--outfile', type = str, help = 'Output file.', required = True)
    args = parser.parse_args()

    distances = getdistances(args.events, args.annotations)
    outfh = open(args.outfile, 'w')
    for event in distances:
        outfh.write(event + '\t' + str(distances[event]) + '\n')
    outfh.close()
