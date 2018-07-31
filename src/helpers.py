import numpy as np
import gzip


def pairwise_diversity(allele_counts, sample_size):
    f = allele_counts / sample_size
    return 2*f*(1-f) / (1 - (1/sample_size))


def folded_sfs(minor_allele_counts, sample_size):
    bins=np.arange(1, sample_size//2 + 2)
    return np.histogram(minor_allele_counts, bins=bins)[0]


def sfs2pi(sfs, sample_size):
    counts = np.arange(1, sample_size//2 + 1)
    pi_weights = pairwise_diversity(counts, sample_size)
    return pi_weights.dot(sfs)


def cross_correlation(X, Y, lim):
    # return np.correlate(X, Y[:-lim], 'valid') / (X.shape[0] - lim)
    n_comps = X.shape[0] - lim
    xy = np.correlate(X, Y[:-lim], 'valid') / n_comps
    x = np.convolve(X, np.ones(n_comps), 'valid') / n_comps
    y = np.nansum(Y[:-lim]) / n_comps
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


def jsfs2corr(msfs, jsfs, n_samples, normalize=True):
    '''Convert marginal and joint site frequency spectra
    to correlations in summary statistics.

    Keyword arguments:
    msfs -- the marginal (average) site frequency spectrum
    jsfs -- the joint (two-site) site frequency spectrum
    n_samples -- the sample size
    normalize -- if True, normalize correlations by product of means
                 if False, return raw correlations
    '''
    n_samples_fold = (n_samples+1)//2
    pi_weights = pairwise_diversity(np.arange(1, n_samples_fold+1), n_samples)

    # Calculate correlation in pi at the two sites
    pi = sfs2pi(msfs, n_samples)
    pi_sq = np.sum(jsfs * pi_weights[:,None] * pi_weights[None,:])
    if normalize:
        pi_corr = (pi_sq/pi**2) - 1
    else:
        pi_corr = pi_sq - pi**2

    lolo_corr = np.zeros(n_samples_fold-1)
    lohi_corr = np.zeros(n_samples_fold-1)
    hihi_corr = np.zeros(n_samples_fold-1)
    # Calculate correlation in reduced jsfs for each cutoff
    for cutoff in np.arange(1, n_samples_fold):
        i = cutoff - 1
        lo = np.sum(msfs[:cutoff])
        hi = np.sum(msfs[cutoff:])
        lolo = np.sum(jsfs[:cutoff, :cutoff])
        lohi = np.sum(jsfs[:cutoff, cutoff:])
        hihi = np.sum(jsfs[cutoff:, cutoff:])
        if normalize:
            lolo_corr[i] = lolo/(lo**2) - 1
            lohi_corr[i] = lohi/(lo*hi) - 1
            hihi_corr[i] = hihi/(hi**2) - 1
        else:
            lolo_corr[i] = lolo - lo**2
            lohi_corr[i] = lohi - lo*hi
            hihi_corr[i] = hihi - hi**2

    return pi_corr, lolo_corr, lohi_corr, hihi_corr


def normalizecorr(readcorrout, n_samples):
    '''Normalize the correlation functions by the product means

    Keyword arguments:
    readcorrout -- tuple output by readcorr
    '''

    sfs, pi_corr, lolo_corr, lohi_corr, hihi_corr = readcorrout
    pi = sfs2pi(sfs, n_samples)
    pi_corr /= pi**2
    for i in range(lolo_corr.shape[0]):
        cutoff = i+1
        lo = np.sum(sfs[:cutoff])
        hi = np.sum(sfs[cutoff:])
        lolo_corr[i] /= lo*lo
        lohi_corr[i] /= lo*hi
        hihi_corr[i] /= hi*hi
    return pi, sfs, pi_corr, lolo_corr, lohi_corr, hihi_corr


### SLiM I/O


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


def readcorr_normed(fn, n_samples):
    '''Read data from a correlation file and normalize correlations'''
    data = readcorr(fn)
    return normalizecorr(data, n_samples)


### msprime I/O


def import_msprime_sfs(file_list, n_samples):
    '''Import pi and the marginal and joint folded sfs from msprime output.

    Keyword arguments:
    file_list -- list of input files containing unfolded
                 marginal and joint sfs
    n_samples -- the sample size of the simulations
    '''

    n_files = len(file_list)
    mSFS = np.zeros((n_files, n_samples-1))
    jSFS_triu = np.zeros((n_files, n_samples*(n_samples-1)//2))

    # Import data from files
    for i, f in enumerate(file_list):
        with open(f) as datafile:
            for line in datafile:
                # Skip header lines
                if line.startswith('#'):
                    continue
                # First two non-header lines contain
                # the marginal and joint SFS.
                mSFS[i,:] = np.array(line.split(), dtype=float)
                jSFS_triu[i,:] = np.array(datafile.readline().split())
                break

    # Unpack the joint SFS from 1D to 2D array
    jSFS = np.zeros((n_files, n_samples-1, n_samples-1))
    for i in range(n_files):
        jSFS[i,:,:][np.triu_indices(n_samples-1)] = jSFS_triu[i,:]
        # Don't double-count the diagonal values
        jSFS[i,:,:][np.diag_indices(n_samples-1)] /= 2
        # Symmetrize distribution
    jSFS += np.transpose(jSFS, axes=(0,2,1))

    # Fold marginal and joint SFS
    mSFS_fold = (mSFS + mSFS[:,::-1])[:,:n_samples//2]
    jSFS_fold = (jSFS + jSFS[:,::-1,:]
                 + jSFS[:,:,::-1]
                 + jSFS[:,::-1,::-1])[:,:n_samples//2,:n_samples//2]
    # Don't double-count the n//2 = n - n//2 values
    if n_samples % 2 == 0:
        mSFS_fold[:,-1] /= 2
        jSFS_fold[:,-1,:-1] /= 2
        jSFS_fold[:,:-1,-1] /= 2
        jSFS_fold[:,-1,-1] /= 4

    pi = sfs2pi(mSFS_fold.T, n_samples)

    return pi, mSFS_fold, jSFS_fold

def import_msprime_corr(file_list, n_samples, normalize=True):
    '''Import pi and the marginal and joint folded sfs from msprime output
    and calculate correlations in pi and the sfs.

    Keyword arguments:
    file_list -- list of input files containing unfolded
                 marginal and joint sfs
    n_samples -- the sample size of the simulations
    normalize -- if True, normalize correlations by product of means
                 if False, return raw correlations
    '''
    n_files = len(file_list)
    n_samples_fold = (n_samples+1)//2

    # Import sfs
    pi, sfs, jsfs = import_msprime_sfs(file_list, n_samples)

    # Calculate correlations
    pi_corr = np.zeros(n_files)
    lolo_corr = np.zeros((n_files, n_samples_fold - 1))
    lohi_corr = np.zeros((n_files, n_samples_fold - 1))
    hihi_corr = np.zeros((n_files, n_samples_fold - 1))
    data = (pi_corr, lolo_corr, lohi_corr, hihi_corr)
    for i in range(n_files):
        new_data = jsfs2corr(sfs[i], jsfs[i], n_samples, normalize=normalize)
        for old, new in zip(data, new_data):
            old[i] = new
    return pi, sfs, jsfs, pi_corr, lolo_corr, lohi_corr, hihi_corr

# fastNeutrino I/O
def get_sfs_from_fastNeutrino(log_fn):
    '''Get expected and observed SFS from fastneutrino log'''
    with open(log_fn) as infile:
        for line in infile:
            if line == 'Expected folded spectrum under inferred demographic model:\n':
                es_line = infile.readline()
            elif line == 'Observed folded spectrum in data:\n':
                os_line = infile.readline()
            elif line.startswith('KL divergence (observed || expected) ='):
                kl_divergence = line.split('=')[1].strip()
    expected_spectrum = np.array(es_line.split(), dtype=float)
    observed_spectrum = np.array(os_line.split(), dtype=float)
    return expected_spectrum, observed_spectrum, kl_divergence
