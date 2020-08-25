import argparse

def fixpseudoT(pseudoTreactivities, outfile):
    reactivitiesfh = open(pseudoTreactivities, 'r')
    outfh = open(outfile, 'w')
    outfh.write(('\t').join(['sequence','rt_start','five_prime_offset','nucleotide','treated_mods','untreated_mods','beta','theta','c']) + '\n')
    currentRNAreacts = []
    for line in reactivitiesfh:
        line = line.strip().split('\t')
        sequence = line[0]
        rt_start = line[1]
        five_prime_offset = line[2]
        nucleotide = line[3]
        treated_mods = line[4]
        untreated_mods = line[5]
        beta = line[6]
        theta = line[7]
        c = line[8]
        if sequence == 'sequence':
            continue
        currentRNAreacts.append(line)

        if int(five_prime_offset) < 131:
            continue
        elif int(five_prime_offset) == 131: #if this is the last nt in the sequence
            nt0_treated_mods = int(currentRNAreacts[0][4])
            nt0_untreated_mods = int(currentRNAreacts[0][5]) #these are the mods at the 0 nt (the *)
            nt1_treated_mods = int(currentRNAreacts[1][4])
            nt1_untreated_mods = int(currentRNAreacts[1][5]) #these are the mods at the 1 nt (the pseduoT)
            new_nt0_treated_mods = str(nt0_treated_mods + nt1_treated_mods) 
            new_nt0_untreated_mods = str(nt0_untreated_mods + nt1_untreated_mods) #new nt0 mods = oldnt0 mods + old nt1 mods
            
            #Write the new nt 0 line.  These always have beta and theta values of -.  
            outfh.write(('\t').join([sequence, '131', '0', '*', new_nt0_treated_mods, new_nt0_untreated_mods, '-', '-', c]) + '\n')
            
            for nt in currentRNAreacts:
                sequence = nt[0]
                rt_start = int(nt[1])
                five_prime_offset = int(nt[2])
                nucleotide = nt[3]
                treated_mods = nt[4]
                untreated_mods = nt[5]
                beta = nt[6]
                theta = nt[7]
                c = nt[8]
                if five_prime_offset == 0 or five_prime_offset == 1: #first line written will be the G at five_prime_offset 2
                    continue
                elif five_prime_offset > 1:
                    outfh.write(('\t').join([sequence, str(rt_start-1), str(five_prime_offset-1), nucleotide, treated_mods, untreated_mods, beta, theta, c]) + '\n')

            currentRNAreacts = [] #reset list

    outfh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pseudoTreactivities', type = str, required = True, help = 'reactivities.out from a spats run where pseudoT mappings were allowed.')
    parser.add_argument('--outfile', type = str, required = True, help = 'Output file of reactivities with pseudoT removed.')
    args = parser.parse_args()

    fixpseudoT(args.pseudoTreactivities, args.outfile)
