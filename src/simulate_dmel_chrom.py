import numpy as np
import msprime
import tskit
import pandas as pd
from sys import argv


def sites(bl_calc, sample_set, site_positions):
    '''
    Computes the expected *derived* (unfolded) site frequency spectrum,
    based on tree lengths, separately in each window.

    :param list sample_set: A list of IDs of samples of length n.
    :param iterable windows: The breakpoints of the windows (including start
        and end, so has one more entry than number of windows).
    :return: A list of lists of length n, one for each window, whose kth
        entry gives the total length of any branches in the marginal trees
        over that window that are ancestral to exactly k of the samples,
        divided by the length of the window.
    '''
    if ((not isinstance(sample_set, list)) or
       len(sample_set) != len(set(sample_set))):
        raise ValueError(
            "elements of sample_sets must be lists without repeated elements.")
    if len(sample_set) == 0:
        raise ValueError("elements of sample_sets cannot be empty.")
    for u in sample_set:
        if not bl_calc.tree_sequence.node(u).is_sample():
            raise ValueError("Not all elements of sample_sets are samples.")

    num_sites = len(site_positions)

    for k in range(num_sites-1):
        if site_positions[k + 1] <= site_positions[k]:
            raise ValueError("Site positions must be increasing.")
    n_out = len(sample_set)
    S = [[0.0 for j in range(n_out)] for _ in range(num_sites)]
    L = [0.0 for j in range(n_out)]
    N = bl_calc.tree_sequence.num_nodes
    X = [int(u in sample_set) for u in range(N)]
    # we will essentially construct the tree
    pi = [-1 for j in range(N)]
    node_time = [bl_calc.tree_sequence.node(u).time for u in range(N)]
    # keep track of where we are for the windows
    chrom_pos = 0.0
    # index of *left-hand* end of the current window
    site_num = 0
    for interval, records_out, records_in in bl_calc.tree_sequence.edge_diffs():
        length = interval[1] - interval[0]
        for sign, records in ((-1, records_out), (+1, records_in)):
            for edge in records:
                dx = 0
                if sign == +1:
                    pi[edge.child] = edge.parent
                dx += sign * X[edge.child]
                dt = (node_time[pi[edge.child]] - node_time[edge.child])
                if X[edge.child] > 0:
                    L[X[edge.child] - 1] += sign * dt
                if sign == -1:
                    pi[edge.child] = -1
                old_X = X[edge.parent]
                X[edge.parent] += dx
                if pi[edge.parent] != -1:
                    dt = (node_time[pi[edge.parent]] - node_time[edge.parent])
                    if X[edge.parent] > 0:
                        L[X[edge.parent] - 1] += dt
                    if old_X > 0:
                        L[old_X - 1] -= dt
                # propagate change up the tree
                u = pi[edge.parent]
                if u != -1:
                    next_u = pi[u]
                    while u != -1:
                        old_X = X[u]
                        X[u] += dx
                        # need to update X for the root,
                        # but the root does not have a branch length
                        if next_u != -1:
                            dt = (node_time[pi[u]] - node_time[u])
                            if X[u] > 0:
                                L[X[u] - 1] += dt
                            if old_X > 0:
                                L[old_X - 1] -= dt
                        u = next_u
                        next_u = pi[next_u]
        while chrom_pos + length >= site_positions[site_num]:
            for j in range(n_out):
                S[site_num][j] += L[j]
            site_num += 1
            if site_num == num_sites:
                return S
        else:
            chrom_pos += length
    return S


# Set simulation parameters
n = 100
output_fn = argv[1]
r = float(argv[2])
model_str = argv[3]
# Optional seed
if len(argv) == 5:
    np.random.seed(int(argv[4]))

# Simulate 10 MB in 1 MB segments
L = int(1e6)
n_reps = 10

# Run simulations
if model_str == 'constant':
    sim = msprime.simulate(model=msprime.StandardCoalescent(),
                           Ne=1/2,
                           recombination_rate=r,
                           length=L,
                           sample_size=n,
                           num_replicates=n_reps)
elif model_str == 'dmel':

    demographic_events = [
                    msprime.PopulationParametersChange(0,
                                                       initial_size=1075793),
                    msprime.PopulationParametersChange(27756,
                                                       initial_size=462688),
                    msprime.PopulationParametersChange(398814,
                                                       initial_size=300000)]
    t2 = 1407233/2
    sim = msprime.simulate(demographic_events=demographic_events,
                           recombination_rate=r/t2,
                           length=L,
                           sample_size=n,
                           num_replicates=n_reps)
else:
    exit(1)

# Get positions of 4D sites from Chr2L
fourD_sites = pd.read_csv('data/dmel-4Dsites.txt.gz', header=None,
                          names=['chr', 'pos'], sep='\t')
# Starting position of central window
start = int(1e6) + 1
positions = np.array(fourD_sites.pos[fourD_sites.chr == '2L'] - start)

branch_lengths = []
for i, rep in enumerate(sim):
    # Get the sites and normalize them
    pos = positions[np.logical_and(positions >= i*L, positions < (i+1)*L)]
    pos -= i*L

    # Compute branch lengths at site positions
    sample_set = list(range(n))
    calculator = tskit.BranchLengthStatCalculator(rep)
    branch_lengths.append(np.array(sites(calculator, sample_set, pos)))

# Save array of branch lengths
np.savez_compressed(output_fn, *branch_lengths)
