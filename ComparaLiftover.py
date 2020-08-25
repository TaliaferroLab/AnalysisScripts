#Gets syntenic regions for every region of a gff.  Currently set up for mouse 3' UTRs and to look in rat, human, cow, and dog.
#Usage: python ComparaLiftover.py <input mouse GFF>

#####Strand info is wrong alot of the time!#######

from cogent.db.ensembl import HostAccount, Genome, Compara
import gffutils
import os
import sys

def ComparaLiftover(gff):
    account = HostAccount('sugarman', 'ensembl', 'ensembl')
    compara = Compara(['mouse', 'rat', 'human', 'cow', 'dog'], Release=61, account=account)
    ratgff = []
    humangff = []
    cowgff = []
    doggff= []

    #Make gff database
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)
    UTRs = db.features_of_type('3\'UTR')

    for UTR in UTRs:
        #Remove stop codons and last 50 nt of UTR
        if UTR.strand == '+':
            UTRstart = int(UTR.start) + 3
            UTRstop = int(UTR.stop) - 50
        elif UTR.strand == '-':
            UTRstart = int(UTR.start) + 50
            UTRstop = int(UTR.stop) - 3
            
        for synt_region in compara.getSyntenicRegions(Species='mouse', CoordName=UTR.chrom.replace('chr',''), Start=UTRstart, End=UTRstop, Strand=UTR.strand, ensembl_coord=True, align_method='PECAN', align_clade='19 amniota vertebrates Pecan'):
            for region in synt_region.Members:
                if region.Region:
                    locdata = str(region.Region.Location).replace('-',':', 1).replace(' ','_').split(':')
                    species = str(locdata[0])
                    chrm = 'chr' + str(locdata[2])
                    start = locdata[3]
                    stop = locdata[4]
                    if str(locdata[5]) == '1' and UTR.strand == '+':
                        strand = '+'
                    elif str(locdata[5]) == '-1' and UTR.strand == '+':
                        strand = '-'
                    elif str(locdata[5]) == '-1' and UTR.strand == '-':
                        strand = '+'
                    elif str(locdata[5]) == '1' and UTR.strand == '-':
                        strand = '-'
                    ID = (str(UTR.id) + '_' + species)
                    if species == 'Rattus_norvegicus':
                        ratgff.append([chrm, 'ALE', '3\'UTR', start, stop, '.', strand, '.', ID])
                    elif species == 'Homo_sapiens':
                        humangff.append([chrm, 'ALE', '3\'UTR', start, stop, '.', strand, '.', ID])
                    elif species == 'Bos_taurus':
                        cowgff.append([chrm, 'ALE', '3\'UTR', start, stop, '.', strand, '.', ID])
                    elif species == 'Canis_familiaris':
                        doggff.append([chrm, 'ALE', '3\'UTR', start, stop, '.', strand, '.', ID])
                

    os.remove(db_fn)
    sys.stderr.write('Succesfully found matches in {0} rat regions, {1} human regions, {2} cow regions and {3} dog regions.\n'.format(len(ratgff), len(humangff), len(cowgff), len(doggff)))
    return ratgff, humangff, cowgff, doggff

if __name__ == '__main__':
    ratfh = open('ALEDistal0.1AxonDISTALUTRs_filtered_rn4.gff3', 'w')
    humanfh = open('ALEDistal0.1AxonDISTALUTRs_filtered_hg19.gff3', 'w')
    cowfh = open('ALEDistal0.1AxonDISTALUTRs_filtered_bosTau4.gff3', 'w')
    dogfh = open('ALEDistal0.1AxonDISTALUTRs_filtered_CanFam2.gff3', 'w')
    ratgff, humangff, cowgff, doggff = ComparaLiftover(sys.argv[1])

    for line in ratgff:
        ratfh.write(('\t').join(line) + '\n')
    for line in humangff:
        humanfh.write(('\t').join(line) + '\n')
    for line in cowgff:
        cowfh.write(('\t').join(line) + '\n')
    for line in doggff:
        dogfh.write(('\t').join(line) + '\n')

    ratfh.close()
    humanfh.close()
    cowfh.close()
    dogfh.close()
    
    
