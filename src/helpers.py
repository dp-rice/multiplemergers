import numpy as np
import gzip

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
    # return np.correlate(X, Y[:-lim], 'valid') / (X.shape[0] - lim)
    n_comps = X.shape[0] - lim
    xy = np.correlate(X, Y[:-lim], 'valid') / n_comps
    x = np.convolve(X, np.ones(n_comps), 'valid') / n_comps
    y = np.sum(Y[:-lim]) / n_comps
    return (xy - x*y)

def windowed_sfs(positions, minor_allele_counts, sample_size, L, window_size):
    sfs_w = np.zeros(((sample_size+1)//2, L//window_size), dtype=int)
    pos_w = positions // window_size
    bins = np.arange(L//window_size + 1)
    for i in range(sfs_w.shape[0]):
        sfs_w[i,:] = np.histogram(pos_w[minor_allele_counts == i+1], bins=bins)[0]
    return sfs_w

def smooth(x, window_len=11, window='hamming'):
    if window_len < 3:
        return x
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat':
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    # Trim off the ends, which are outside the x-range of the data
    trim = window_len//2
    return y[trim:-trim]

# This is faster than numpy loadtxt because it pre-allocates the space
def loadints(fn, rows, cols): 
    ret = np.zeros((rows, cols), dtype=int)
    with gzip.open(fn, 'rb') as infile:
        i = 0
        for line in infile:
            if line.startswith(b'#'):
                continue
            ret[i,:] = line.split()
            i += 1
    return ret

def readcorr(fn):
    with gzip.open(fn, 'rb') as infile:
        # SFS header and SFS
        infile.readline()
        sfs = np.array([float(x) for x in infile.readline().split()])

        # PI_CORR header and PI_CORR
        infile.readline()
        pi_corr = np.array([float(x) for x in infile.readline().split()])

        lolo_corr = np.zeros((len(sfs)-1, len(pi_corr)))
        lohi_corr = np.zeros((len(sfs)-1, len(pi_corr)))
        hihi_corr = np.zeros((len(sfs)-1, len(pi_corr)))

        # LOLO_CORR header and LOLO_CORR
        infile.readline()
        for i in range(lolo_corr.shape[0]):
            line = infile.readline()
            y = np.array([float(x) for x in line.split()])
            lolo_corr[i,:] = y
        # LOHI_CORR header and LOHI_CORR
        infile.readline()
        for i in range(lohi_corr.shape[0]):
            line = infile.readline()
            y = np.array([float(x) for x in line.split()])
            lohi_corr[i,:] = y
        # HIHI_CORR header and HIHI_CORR
        infile.readline()
        for i in range(hihi_corr.shape[0]):
            line = infile.readline()
            y = np.array([float(x) for x in line.split()])
            hihi_corr[i,:] = y

        return sfs, pi_corr, lolo_corr, lohi_corr, hihi_corr
