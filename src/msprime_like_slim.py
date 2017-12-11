import sys
import numpy as np
import argparse

#FIXME: this is a hack to use my local version of msprime
# sys.path.insert(1, '/users/danielrice/msprime-lambda/')
# sys.path.insert(1, '/home/dpr/mmc_genomics/src/msprime-lambda/')
sys.path.insert(1, '/Users/dpr/mmc_genomics/src/msprime-lambda/')
import msprime
# import jsfs

# FIXME: take these from a config file
n_samples = 100
N = 500
L = int(1e8)
mu_n = 1e-7
r = 1e-8

# Run simulation
treeseq = msprime.simulate(
        sample_size=n_samples,
        length=L,
        recombination_rate=r,
        Ne=N,
        mutation_rate=mu_n)

# Strings to format the output like SLiM
header_str = '#OUT: {} . . {}\nMutations:\n'
line_str = '. . m0 {} . . . . {}\n'
sys.stdout.write(header_str.format(N, n_samples))
# Loop over trees in the tree sequence, getting the positions and allele counts
for tree in treeseq.trees():
    for mutation in tree.mutations():
        position = int(np.floor(mutation.position))
        allele_count = tree.get_num_leaves(mutation.node)
        sys.stdout.write(line_str.format(position, allele_count))
# This is so the parsing script knows to stop
sys.stdout.write('Genomes:\n')
