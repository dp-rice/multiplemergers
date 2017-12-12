import sys
import numpy as np
import helpers as h

infn = sys.argv[1]
L_str = sys.argv[2]
w_str = sys.argv[3]
outfn = sys.argv[4]

if 'e' in L_str:
    L = int(float(L_str))
else:
    L = int(L_str)

if 'e' in w_str:
    window_size = int(float(w_str))
else:
    window_size = int(w_str)

# sys.stderr.write('Importing data...\n')
sample_size, positions, allele_counts = h.import_slim_output(infn)
minor_allele_counts = np.minimum(allele_counts, sample_size - allele_counts)

# sys.stderr.write('Computing windowed statistics...\n')
sfs_w = h.windowed_sfs(positions, minor_allele_counts, sample_size, L, window_size)

# sys.stderr.write('Saving to file...\n')
np.savetxt(outfn, sfs_w, fmt='%d', header='window_size={}'.format(window_size))
