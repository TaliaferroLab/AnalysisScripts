import sys

def MISOtoTable(misofile, eventtype, outfile):
    infh = open(misofile, 'r')
    lines = []
    for line in infh:
        line = line.strip().split('\t')
        if line[0] == 'event_name':
            continue
        elif line[0] != 'event_name':
            eventname = line[0]
            sample1psi = line[1]
            sample2psi = line[4]
            deltapsi = line[7]
            lines.append([eventname, sample1psi, sample2psi, deltapsi, eventtype])
    infh.close()

    outfh = open(outfile, 'w')
    for line in lines:
        outfh.write(('\t').join(line) + '\n')

    outfh.close()

if __name__ == '__main__':
    MISOtoTable(sys.argv[1], sys.argv[2], sys.argv[3])
