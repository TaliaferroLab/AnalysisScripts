import sys
import os

def formisosummarycluster(directoryofsummaries):
    for file in os.listdir(directoryofsummaries):
        infh = open(file, 'r')
        outfh = open(file[:-3] + 'forcluster.bed', 'w')
        for line in infh:
            line = line.strip().split('\t')
            psi = float(line[4])
            confint = line[5]
            if confint != 'NA':
                '''
                if float(confint) <= 0.3: #confidence interval cutoff
                    outfh.write(line[0] + ':' + line[1] + ':' + line[2] + ':' + line[3] + '\t' + line[4] + '\n')
                '''
                if float(confint) <= 0.3:
                    outfh.write(line[0] + ':' + line[1] + ':' + line[2] + ':' + line[3] + '\t' + line[5] + '\n')
        infh.close()
        outfh.close()

formisosummarycluster(sys.argv[1])
