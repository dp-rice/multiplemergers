import sys
import gzip
import numpy as np

n_samples = int(sys.argv[1])
# The total length of the simulated chromosome
L_str = sys.argv[2]
# The window size
w_str = sys.argv[3]
# The remaining args are a list of filenames
filenames = sys.argv[4:]

if 'e' in L_str:
    L = int(float(L_str))
else:
    L = int(L_str)

if 'e' in w_str:
    w = int(float(w_str))
else:
    w = int(w_str)

lim = L//(10*w)
n_files = len(filenames)

sfs = np.zeros((n_samples+1)//2)
pi_corr = np.zeros(lim + 1)
lohi_corr = np.zeros((len(sfs)-1, len(pi_corr)))

for fn in filenames:
    with gzip.open(fn, 'rb') as infile:
        # SFS header and SFS
        infile.readline()
        sfs += np.array([float(x) for x in infile.readline().split()])
        # sfs += np.array(infile.readline().split())
        # PI_CORR header and PI_CORR
        infile.readline()
        pi_corr += np.array([float(x) for x in infile.readline().split()])
        # LOHI_CORR header and LOHI_CORR
        infile.readline()
        for i, line in enumerate(infile):
            y = np.array([float(x) for x in line.split()])
            lohi_corr[i,:] += y

sfs /= n_files
pi_corr /= n_files 
lohi_corr /= n_files

# Write output.
sys.stdout.write('#SFS:\n')
sys.stdout.write(' '.join(map(str, sfs)) + '\n')

sys.stdout.write('#PI_CORR:\n')
sys.stdout.write(' '.join(map(str, pi_corr)) + '\n')

sys.stdout.write('#LOHI_CORR:\n')
for i in range(lohi_corr.shape[0]):
    sys.stdout.write(' '.join(map(str, lohi_corr[i,:])) + '\n')
