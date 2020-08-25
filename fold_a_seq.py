import subprocess
import sys

def fold_a_seq(seq):
    """
    Returns the folding energy for a sequence
    THis is much slower than doing it in a batch.
    """
    command = 'RNAfold'
    job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    job.stdin.write(seq)
    output = job.communicate()
    return float(output[0].replace('(', ' ').replace(')', ' ').split()[-1])


if __name__ == '__main__':
    print fold_a_seq(sys.argv[1])
