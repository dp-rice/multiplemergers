import argparse
import numpy as np
import helpers as h
from scipy.stats import hypergeom

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
    # If fewer than n_downsample observations, return zeros.
    if n_downsample > n_obs:
        sfs = np.zeros(n_downsample+1)
    # Otherwise, use hypergeometric probability as expected sfs
    else:
        x = np.arange(0, n_downsample+1)
        sfs = hypergeom(n_obs, allele_count, n_downsample).pmf(x)
    return sfs

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
w = args.window_size
n_ds = args.n_downsample

# Strip off the "chr" from chromosme names.
chrom = args.chromosome[3:]
# Import 4-fold degenerate sites
with open(args.fourfold_file) as infile:
    # Convert 4d sites from one-indexed to zero-indexed.
    fourDsites = [int(line.split()[1]) - 1 for line in infile
                    if line.startswith(chrom)]

# Import allele-count data
data = h.loadints(args.data_file, args.chrom_len, 2)

# Initialize per-locus SFS matrix
sfs = np.zeros((n_ds+1, args.chrom_len // w + 1))
# Calculate downsampled sfs
for pos in fourDsites:
    sfs[:, pos//w] +=  downsample_sfs(data[pos,0], data[pos,1], n_ds)

# Print header
print(n_ds, args.chrom_len // w + 1)
# Print SFS
for row in sfs:
    print(' '.join(map(str,row)))
