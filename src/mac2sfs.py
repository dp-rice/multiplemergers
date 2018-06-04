#import sys
import argparse
import numpy as np
import helpers as h

def downsample_sfs(n_obs, allele_count, n_downsample):
    """
    Calculate expected site frequency spectrum for a smaller sample size.

    Arguments:
    n_obs -- the original number of haploid genotypes
    minor_allele_count -- the original allele count
    n_downsample -- the new (smaller) sample size

    Returns:
    sfs -- the expected downsampled SFS as a length--n_downsample + 1 numpy array
    """
    pass

parser = argparse.ArgumentParser(description='''Calculate windowed
                    and downsampled site frequency spectrum.''')
parser.add_argument('chromosome', metavar='c', type=str,
                    help='Chromosome name')
parser.add_argument('chrom_len', metavar='L', type=int,
                    help='Chromsome length')
parser.add_argument('window_size', metavar='w', type=int,
                    help='Length of windows for binned SFS')
parser.add_argument('n_downsample', metavar='nd', type=int,
                    help='New sample size')
parser.add_argument('data_file', type=str,
                    help='File containing minor allele count at each site')
parser.add_argument('fourfold_file', type=str,
                    help='File containing a list of four-fold degenerate sites')
args = parser.parse_args()

# Strip off the "chr" from chromosme names.
chrom = args.chromosome[3:]

# Import 4-fold degenerate sites
with open(args.fourfold_file) as infile:
    # Convert 4d sites from one-indexed to zero-indexed.
    fourDsites = [int(line.split()[1]) - 1 for line in infile
                    if line.startswith(chrom)]

# Import allele-count data
data = h.loadints(args.data_file, args.chrom_len, 2)

print(fourDsites[:10])
print(data.shape)

# Initialize per-locus SFS matrix

# Calculate downsampled sfs

# Print header

# Print SFS
