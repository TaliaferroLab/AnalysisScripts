import sys

def trimto28(fastq, outfastq):
    lineindexer = 1
    counter = 0
    fastqfh = open(fastq, 'r')
    outfh = open(outfastq, 'w')
    outfh.close()
    for line in fastqfh:
        if lineindexer == 1:
            with open(outfastq, 'a') as outfh:
                outfh.write(line)
            lineindexer +=1
            continue
        elif lineindexer == 2:
            with open(outfastq, 'a') as outfh:
                outfh.write(line[:28] + '\n')
            lineindexer +=1
            continue
        elif lineindexer == 3:
            with open(outfastq, 'a') as outfh:
                outfh.write(line)
            lineindexer +=1
            continue
        elif lineindexer == 4:
            with open(outfastq, 'a') as outfh:
                outfh.write(line[:28] + '\n')
            lineindexer = 1
            counter +=1
            if counter %10000000 == 0:
                print 'Considering read {0}'.format(counter)
            continue

trimto28(sys.argv[1], sys.argv[2])
