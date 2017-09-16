import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("nSamples", type=int, help="number of samples to subsample to. Must be less than or equal to the number of samples in the file")
args = parser.parse_args()

alleles = np.zeros(args.nSamples, dtype='S1')
for i_line, line in enumerate(sys.stdin):
    data = list(line.strip()[:args.nSamples])
    # Masked site or site with uncalled genotypes
    if 'N' in data:
        maf = 'NaN'
    # Homozygous site
    elif len(data) == 1:
        maf = '0'
    else:
        alleles[:] = data
        bases, counts = np.unique(alleles, return_counts=True)
        # Multi-allelic site
        if len(counts) > 2:
            maf = 'NaN'
        # Polymorphic site
        else:
            maf = str(np.min(counts)%args.nSamples)
    sys.stdout.write(maf + '\n')
