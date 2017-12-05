import numpy as np

def pairwise_diversity(allele_counts, sample_size):
    f = allele_counts / sample_size
    return 2*f*(1-f) / (1 - (1/sample_size))

def folded_sfs(minor_allele_counts, sample_size):
    bins=np.arange(1, (sample_size+1)//2 + 2)
    return np.histogram(minor_allele_counts, bins=bins)[0]

def sfs2pi(sfs, sample_size):
    counts = np.arange(1, (sample_size+1)//2 + 1)
    pi_weights = pairwise_diversity(counts, sample_size)
    return pi_weights.dot(sfs)

def import_slim_output(filename):
    positions = []
    allele_counts = []
    with open(filename) as infile:
        line = ''
        # Scan to the begining of the sample output
        while not line.startswith('#OUT'):
            line = infile.readline()
        sample_size = int(line.split()[-1])
        while not line.startswith('Mutations:'):
            line = infile.readline()
        # Get the segregating mutations and allele counts
        while infile:
            line = infile.readline()
            if line.startswith('Genomes:'):
                break
            sline = line.split()
            ac = int(sline[8])
            if ac < sample_size:
                positions.append(int(sline[3]))
                allele_counts.append(ac)
    return sample_size, np.array(positions), np.array(allele_counts)

def cross_correlation(X, Y, lim):
    return np.correlate(X, Y[:-lim], 'valid') / (X.shape[0] - lim)

def windowed_sfs(positions, minor_allele_counts, sample_size, L, window_size):
    sfs_w = np.zeros(((sample_size+1)//2, L//window_size), dtype=int)
    pos_w = positions // window_size
    bins = np.arange(L//window_size + 1)
    for i in range(sfs_w.shape[0]):
        sfs_w[i,:] = np.histogram(pos_w[minor_allele_counts == i+1], bins=bins)[0]
    return sfs_w
