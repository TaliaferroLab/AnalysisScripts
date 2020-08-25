#Takes a list of MISO Events and a MISO annotation file.  Retrieves MISO annotations that are in the events file.  Sorting of annotations is lost.  May need to re-sort ouputted gff3 using ReorderGff3.py.

#Usage python RetrieveMISOAnnotations.py <events> <annotation> <outfile>

import gffutils
import os
import sys

annotation_fn = sys.argv[2]
db_fn = os.path.basename(annotation_fn) + '.db'
outfh = open(sys.argv[3],'w')

if os.path.isfile(db_fn) == False:
    gffutils.create_db(annotation_fn, db_fn)

db = gffutils.FeatureDB(db_fn)

eventsfh = open(sys.argv[1],'r')

def getEventIntervals(annotation_fn, eventsfh):

    events = []
    for line in eventsfh:
        line = line.strip()
        line = line.split('\t')
        events.append(line[0])
       

    genes = db.features_of_type('gene')

    for gene in genes:
       for event in events:
           if event == gene.id:
               outfh.write(str(gene) + '\n')
               for child in db.children(gene):
                   outfh.write(str(child) + '\n')

    outfh.close()



getEventIntervals(annotation_fn, eventsfh)
os.remove(db_fn)
if os.path.isfile(db_fn + '-shm'):
    os.remove(db_fn + '-shm')
if os.path.isfile(db_fn + '-wal'):
    os.remove(db_fn + '-wal')
