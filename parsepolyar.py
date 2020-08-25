import sys
import argparse

def parsepolyar(polyaroutput):
    infh = open(polyaroutput, 'r')
    queryscores = {}
    seqswithpolyAsites = 0
    queries = 0
    for line in infh:
        line = line.strip()
        if 'Query Name' in line:
            queries +=1
            query = line.split('>')[2]
            currentqueryscores = []
            polyAsites = 0
        if 'W:' in line: #line with scores
            polyAsites +=1
            score = float(line.split('W:')[1])
            currentqueryscores.append(score)
        if 'Totally found' in line: #last line of query group
            seqswithpolyAsites +=1
            maxscore = max(currentqueryscores)
            queryscores[query] = [maxscore, polyAsites]
    infh.close()
    sys.stderr.write('PolyAsites were found in {0} of {1} sequences.\n'.format(seqswithpolyAsites, queries))
    return queryscores

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type = str, help = 'Output of polyar to parse.')
    parser.add_argument('--output', type = str, help = 'File in which to output results.')
    args = parser.parse_args()
    queryscores = parsepolyar(args.input)
    outfh = open(args.output, 'w')
    outfh.write('Sequence' + '\t' + 'Max_polyAsite_score' + '\t' + 'Number_of_sites' + '\n')
    for query in queryscores:
        score = queryscores[query][0]
        polyAsites = queryscores[query][1]
        outfh.write(query + '\t' + str(score) + '\t' + str(polyAsites) + '\n')
    outfh.close()
