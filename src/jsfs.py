import sys
import numpy as np
import msprime
from collections import deque

def interval_overlap(inter1, inter2):
    left = max(inter1[0], inter2[0])
    right = min(inter1[1], inter2[1])
    # FIXME: tolerance for rounding?
    if left >= right:
        return None
    else:
        return (left, right)

def interval_length(interval):
    return interval[1] - interval[0]

def merge(interval_lists):
    merged_intervals = interval_lists[0]
    for merging_intervals in interval_lists[1:]:
        new_intervals = []
        for i1, s1 in merged_intervals:
            for i2, s2 in merging_intervals:
                io = interval_overlap(i1,i2)
                if io:
                    new_intervals.append((io, s1+s2))
        merged_intervals = new_intervals
    return merged_intervals

def get_sfs_from_mutations(treeSequence, error_rate):
    sample_size = treeSequence.get_sample_size()
    T = np.zeros(sample_size - 1)
    for tree in treeSequence.trees():
        for mutation in tree.mutations():
            # node = mutation.node
            node = mutation[1]
            num_leaves = tree.get_num_leaves(node)
            if error_rate > 0:
                num_leaves += np.random.poisson((sample_size - num_leaves) * error_rate) - np.random.poisson(num_leaves*error_rate)
                # Throw out sites that would not be identified thanks to errors
                if num_leaves <= 0 or num_leaves >= sample_size:
                    continue
            T[num_leaves - 1] += 1
    return T

def get_sfs(treeSequence):
    # This algorithm loops over coalescent records instead of trees.
    # For each node, it keeps list of the intervals in which the node exists and the number of descendents in each interval
    # Since the records are time-sorted, you can recursively calculate these intervals for each node, moving up the tree

    sample_size = treeSequence.get_sample_size()
    # Set up leaves: time = 0, interval = (0.0,1.0), size = 1
    nodes = [(0,[((0.0, 1.0), 1)])] * sample_size
    # Initialize SFS times
    T = np.zeros(sample_size - 1)
    records = treeSequence.records()
    for rec in records:
        rec_interval = (rec.left, rec.right)
        # For each child, find the get the list of intervals that overlap with this record
        subintervals = []
        for child in rec.children:
            child_time, child_intervals = nodes[child]
            # Branch length
            t = rec.time - child_time
            child_subintervals = []
            for interval, size in child_intervals:
                 io = interval_overlap(rec_interval, interval)
                 if io:
                     # Update T
                     T[size-1] += interval_length(io) * t
                     child_subintervals.append((io, size))
            subintervals.append(child_subintervals)
        # Combine the intervals from each child to get the list of intervals + sizes for the parent node
        merged_subintervals = merge(subintervals)
        # Update nodes
        try:
            node = nodes[rec.node]
        except IndexError:
            node = (rec.time, [])
            nodes.append(node)
        node[1].extend(merged_subintervals)
    return T

def sfs_from_tree(tree, n_samples):
    sfs = np.zeros(n_samples - 1)
    root = tree.get_root()
    for u in tree.nodes():
        if u == root:
            continue
        t = tree.get_branch_length(u)
        l = tree.get_num_leaves(u)
        sfs[l-1] += t
    return sfs

def get_jsfs(tree_sequence):
    n_samples = tree_sequence.get_sample_size()
    # tree = next(tree_sequence.trees())
    # tau1 = sfs_from_tree(tree, n_samples)
    # deq = deque(tree_sequence.trees(), maxlen=1)
    # if len(deq) == 0:
    #     tau2 = np.copy(tau1)
    # else:
    #     tree = deq.pop()
    #     tau2 = sfs_from_tree(tree, n_samples)
    # return tau1, tau2

    n_trees = tree_sequence.get_num_trees()
    for i, tree in enumerate(tree_sequence.trees()):
        if i == 0:
            tau1 = sfs_from_tree(tree, n_samples)
        # NOTE: this could be the same tree if n_trees == 1.
        if i == n_trees - 1:
            tau2 = sfs_from_tree(tree, n_samples)
    return tau1, tau2
    TAU = np.zeros((2, tree_sequence.get_sample_size() - 1))
    for tree in tree_sequence.trees():
        interval_length = tree.get_length()
    return TAU

def get_sfs_naive(tree_sequence):
    TAU = np.zeros(tree_sequence.get_sample_size() - 1)
    for tree in tree_sequence.trees():
        interval_length = tree.get_length()
        root = tree.get_root()
        for u in tree.nodes():
            if u == root:
                continue
            t = tree.get_branch_length(u)
            l = tree.get_num_leaves(u)
            # Adjust the branch lengths by the interval lengths to account for mutation opportunity
            TAU[l-1] += t * interval_length
    return TAU


if __name__ == "__main__":
    from time import time

    nLoci = 1000
    nSamples = 100
    r = 1.0
    alpha = 1.5
    TAU = np.zeros((nLoci, nSamples - 1))
    TAU_naive = np.zeros((nLoci, nSamples - 1))

    population_configurations = [msprime.PopulationConfiguration(sample_size=nSamples, multiple_merger_para=alpha)]

    simulations = msprime.simulate(num_replicates=nLoci, recombination_rate=r,
                                    population_configurations=population_configurations)
    t_test = 0
    t_naive = 0
    for i_rep, rep in enumerate(simulations):
        tstart = time()
        TAU[i_rep,:] = get_sfs(rep)
        t_test += time() - tstart

        tstart = time()
        TAU_naive[i_rep,:] = get_sfs_naive(rep)
        t_naive += time() - tstart

    # print 'T_naive:', t_naive
    # print 'T_test: ', t_test
    # print 'Max diff:', np.max(abs(TAU - TAU_naive))
