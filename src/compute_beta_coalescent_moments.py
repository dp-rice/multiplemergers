import numpy as np
import coalescentmoments as cm
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("nSamples", type=int, help="number of samples to simulate")
parser.add_argument("alpha", type=float, help='beta coalescent parameter alpha')
parser.add_argument('--folded_sfs', dest='folded', action='store_true',
        help='Record folded site frequency spectrum (DEFAULT=True)')
parser.add_argument('--unfolded_sfs', dest='folded', action='store_false',
        help='Record unfolded site frequency spectrum (DEFAULT=False)')
parser.set_defaults(folded=True)
args = parser.parse_args()


M1, M2 = cm.sfs_moments(args.nSamples, alpha=args.alpha)

# Print in a format compatible with other scripts
print('#N_SAMPLES={}'.format(args.nSamples))
print('#ALPHA={}'.format(args.alpha))
print('#SFS_FOLDED={}'.format(args.folded))
if args.folded:
    M1_folded, M2_folded = cm.fold_sfs_moments(args.nSamples, M1, M2)
    print(' '.join([str(x) for x in M1_folded]))
    # Print second moment matrix in condensed form (only the upper triangular portion)
    print(' '.join([str(x) for x in M2_folded[np.triu_indices(M2_folded.shape[0])]]))
else:
    print(' '.join([str(x) for x in M1]))
    print(' '.join([str(x) for x in M2[np.triu_indices(M2.shape[0])]]))
