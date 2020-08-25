#Divides the UTR into 100 approximately equally sized chunks (rounding causes some deviations).  Calculates how many clusters hit each chunk in each UTR.  Then calculates densities, which is the number of hits in each chunk, divided by how many times that chunk was interrogated (which is the same as the number of UTRs searched).

#For densities of only one UTR: python 3UTRMetagene.py <gfffile> <clusters bed file> <output.txt>
#For relative densities of two UTRs: python 3UTRMetagene.py <gfffile1> <gfffile2> <clusters bed file> <output.txt>
#Relative densities are gff1 / gff2

import sys
import os
import tabix
import pysam
import gffutils


def threeUTRmetagene(gff, clusters):
    gff_fn = gff
    db_fn = os.path.basename(gff_fn) + '.db'

    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn)

    db = gffutils.FeatureDB(db_fn)

    clustertabix = pysam.Tabixfile(clusters)
    number_of_UTRs = float(len(list(db.features_of_type('3\'UTR'))))
    three_prime_UTRs = db.features_of_type('3\'UTR')
    binhits = {}
    bindensities = {}

    for UTR in three_prime_UTRs:
        binnumber = 1
        lowerpercent = 0
        upperpercent = 1
        UTRlength = float((UTR.stop - UTR.start))
        
        while upperpercent <= 100:
            nt_window = [int(round((UTRlength/100)*lowerpercent)), int(round((UTRlength/100)*upperpercent))]
            windowstart = nt_window[0]
            windowstop = nt_window[1]
            #For every time a cluster bed entry overlaps this particular sequence window
            for bed in clustertabix.fetch(str(UTR.chrom), UTR.start + windowstart, UTR.start + windowstop, parser = pysam.asBed()):
                if UTR.strand == bed.strand:
                    if binhits.has_key(binnumber):
                        binhits[binnumber] +=1
                    else:
                        binhits[binnumber] = 1
            lowerpercent +=1
            upperpercent +=1
            binnumber +=1

    os.remove(db_fn)

    #Number of hits in every unpopulated bin is 0
    for x in range(1,101):
        if x not in binhits:
            binhits[x] = 0
    
    for binnumber in binhits:
        rawhits = binhits[binnumber]
        density = float(rawhits / number_of_UTRs)
        bindensities[binnumber] = density

    return bindensities


#For comparing two UTR gff files 
def threeUTRmetagene_relative_densities(gff1, gff2, clusters):
    gff1_fn = gff1
    gff2_fn = gff2
    db1_fn = os.path.basename(gff1_fn) + '.db'
    db2_fn = os.path.basename(gff2_fn) + '.db'
    
    if os.path.isfile(db1_fn) == False:
        gffutils.create_db(gff1_fn, db1_fn)

    if os.path.isfile(db2_fn) == False:
        gffutils.create_db(gff2_fn, db2_fn)

    db1 = gffutils.FeatureDB(db1_fn)
    db2 = gffutils.FeatureDB(db2_fn)

    clustertabix = pysam.Tabixfile(clusters)
    number_of_UTRs_1 = float(len(list(db1.features_of_type('3\'UTR'))))
    number_of_UTRs_2 = float(len(list(db2.features_of_type('3\'UTR'))))
    three_prime_UTRs_1 = db1.features_of_type('3\'UTR')
    three_prime_UTRs_2 = db2.features_of_type('3\'UTR')
    binhits_1 = {}
    bindensities_1 = {}
    binhits_2 = {}
    bindensities_2 = {}
    relativedensities = {}

    for UTR in three_prime_UTRs_1:
        binnumber = 1
        lowerpercent = 0
        upperpercent = 1
        UTRlength = float((UTR.stop - UTR.start))
        
        while upperpercent <= 100:
            nt_window = [int(round((UTRlength/100)*lowerpercent)), int(round((UTRlength/100)*upperpercent))]
            windowstart = nt_window[0]
            windowstop = nt_window[1]
            for bed in clustertabix.fetch(str(UTR.chrom), UTR.start + windowstart, UTR.start + windowstop, parser = pysam.asBed()):
                if UTR.strand == bed.strand:
                    if binhits_1.has_key(binnumber):
                        binhits_1[binnumber] +=1
                    else:
                        binhits_1[binnumber] = 1
            lowerpercent +=1
            upperpercent +=1
            binnumber +=1

    for UTR in three_prime_UTRs_2:
        binnumber = 1
        lowerpercent = 0
        upperpercent = 1
        UTRlength = float((UTR.stop - UTR.start))
        
        while upperpercent <= 100:
            nt_window = [int(round((UTRlength/100)*lowerpercent)), int(round((UTRlength/100)*upperpercent))]
            windowstart = nt_window[0]
            windowstop = nt_window[1]
            for bed in clustertabix.fetch(str(UTR.chrom), UTR.start + windowstart, UTR.start + windowstop, parser = pysam.asBed()):
                if UTR.strand == bed.strand:
                    if binhits_2.has_key(binnumber):
                        binhits_2[binnumber] +=1
                    else:
                        binhits_2[binnumber] = 1
            lowerpercent +=1
            upperpercent +=1
            binnumber +=1

    os.remove(db1_fn)
    os.remove(db2_fn)

    #Any empty bins should have a density of 0
    for x in range(1,101):
        if x not in binhits_1:
            binhits_1[x] = 0
        if x not in binhits_2:
            binhits_2[x] = 0
            
    for binnumber in binhits_1:
        rawhits = binhits_1[binnumber]
        density = float(rawhits / number_of_UTRs_1)
        bindensities_1[binnumber] = density

    for binnumber in binhits_2:
        rawhits = binhits_2[binnumber]
        density = float(rawhits / number_of_UTRs_2)
        bindensities_2[binnumber] = density

    for binnumber in bindensities_1:
        if bindensities_2[binnumber] != 0:
            relativedensities[binnumber] = float((float(bindensities_1[binnumber]) / float(bindensities_2[binnumber])))
    for binnumber in bindensities_1:
        if bindensities_2[binnumber] == 0:
            relativedensities[binnumber] = relativedensities[min(relativedensities, key=relativedensities.get)] ####FIX THIS###

    return relativedensities
                            
    
if __name__ == '__main__':
    if len(sys.argv) == 4:
        outfh = open(sys.argv[3], 'w')
        outfh.write('bin' + '\t' + 'density' + '\n')
        dic = threeUTRmetagene(sys.argv[1], sys.argv[2])
        for binnumber in dic:
            outfh.write(str(binnumber) + '\t' + str(dic[binnumber]) + '\n')
        outfh.close()

    elif len(sys.argv) == 5:
        outfh = open(sys.argv[4], 'w')
        outfh.write('bin' + '\t' + 'relative_density' + '\n')
        dic = threeUTRmetagene_relative_densities(sys.argv[1], sys.argv[2], sys.argv[3])
        for binnumber in dic:
            outfh.write(str(binnumber) + '\t' + str(dic[binnumber]) + '\n')
        outfh.close()
