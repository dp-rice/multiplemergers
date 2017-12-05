import sys
import numpy as np
import helpers as h

filename = sys.argv[1]
n_samples = int(sys.argv[2])
# The total length of the simulated chromosome
L_str = sys.argv[3]
# The window size
w_str = sys.argv[4]

if 'e' in L_str:
    L = int(float(L_str))
else:
    L = int(L_str)

if 'e' in w_str:
    w = int(float(w_str))
else:
    w = int(w_str)

# The largest distance at which to calculate the cross-correlations
lim = L//(10*w)
# The largest frequency cutoff between low and high frequency mutations
lo_max = (n_samples+1)//2 - 1

# Load windowed site frequency spectrum and calculate average sfs.
wsfs = h.loadints(filename, (n_samples+1)//2, L//w)
sfs = np.sum(wsfs, axis=1) / L

# Calculate windowed pi and its autocorrelation function.
pi_w = h.sfs2pi(wsfs, n_samples) / w
pi_corr = h.cross_correlation(pi_w, pi_w, lim)

# Calculate the cross correlation between high- and low-frequency mutations,
# varying the cutoff between low and high.
lohi_corr = np.zeros((lo_max, lim + 1))
for freq_cutoff in range(1, lo_max + 1):
    lo = np.sum(wsfs[:freq_cutoff,:], axis=0) / w
    hi = np.sum(wsfs[freq_cutoff:,:], axis=0) / w
    lohi_corr[freq_cutoff-1] = h.cross_correlation(lo, hi, lim)

# Write output.
sys.stdout.write('#SFS:\n')
sys.stdout.write(' '.join(map(str, sfs)) + '\n')

sys.stdout.write('#PI_CORR:\n')
sys.stdout.write(' '.join(map(str, pi_corr)) + '\n')

sys.stdout.write('#LOHI_CORR:\n')
for i in range(lo_max):
    sys.stdout.write(' '.join(map(str, lohi_corr[i,:])) + '\n')
