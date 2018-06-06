import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Write fastNeutrino input files with constant-N epochs")
parser.add_argument('ancestral_size', metavar='N_anc', type=float,
                    help='Ancestral population size')
parser.add_argument('num_epochs', metavar='n', type=int,
                    help='Number of constant-size epochs (default=1)')
parser.add_argument('--spacing', type=str, default='free',
                    help='Spacing of the epoch boundaries')

args = parser.parse_args()

# One epoch is the ancestral state.
n = args.num_epochs-1

# TODO: Make these flexible
start = 10
end = 10*args.ancestral_size

if args.spacing == 'free':
    times = ['?t{}'.format(i) for i in range(n)]
elif args.spacing == 'fixed':
    times = np.logspace(np.log10(start), np.log10(end), n)

# First line of the input file is the ancestral size.
print(args.ancestral_size)
for i, t in enumerate(times):
    print('c', '?n{}'.format(i), t)
# There aren't any constrained parameters.
print(0)
print(0)
