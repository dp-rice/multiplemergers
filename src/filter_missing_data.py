import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("nSamples", type=int, help="number of samples in file")
parser.add_argument("maxN", type=int, help='Maximum number of non-genotyped samples allowed to include a site')
parser.add_argument("-p", "--positions_in_file", action="store_true", help="input file contains position numbers")
parser.add_argument("-a", "--maxAlleles", type=int, default=2, help="maximum number of alleles allowed (default=2)")
args = parser.parse_args()

alleles = np.zeros(args.nSamples, dtype='S1')
for i_line, line in enumerate(sys.stdin):
    if args.positions_in_file:
        data = list(line.split()[1])
        outline = line
    else:
        data = list(line.strip())
        outline = str(i_line+1) + '\t' + line

    # Skip masked sites and sites with more than maxN ungenotyped samples
    if data == ['N'] or data.count('N') > args.maxN:
        continue
    # Always write monomorphic sites
    elif len(data) == 1:
        sys.stdout.write(outline)
    # Make sure there are no more than maxAlleles alleles at the site
    else:
        alleles[:] = data
        n_alleles = np.unique(alleles[alleles != 'N']).size
        if n_alleles <= args.maxAlleles:
            sys.stdout.write(outline)
