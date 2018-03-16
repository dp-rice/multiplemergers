import sys
import argparse
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("nSamples", type=int, help="number of samples to subsample to. Must be less than or equal to the number of samples in the file")
args = parser.parse_args()

for i_line, line in enumerate(sys.stdin):
    alleles = line.strip()[:args.nSamples]
    c = Counter(alleles)
    nN = c['N']
    # Remove the Ns so that we can count the number of alleles
    del c['N']
    # Masked site
    if len(c) == 0:
        nob = 0
        mac = 0
    # Monomorphic site
    elif len(c) == 1:
        nob = args.nSamples - nN
        mac = 0
    elif len(c) == 2:
        nob = args.nSamples - nN
        mac = min(c.values())
    # Mask tri-allelic sites
    else:
        nob = 0
        mac = 0
    sys.stdout.write("{} {}\n".format(nob, mac))
