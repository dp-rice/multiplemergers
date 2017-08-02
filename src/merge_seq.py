import sys
import tarfile

# We will read the files in chunks of size CHUNKSIZE bytes
CHUNKSIZE=4096 

#TODO: argparse
tar_fn = sys.argv[1]
sample_list = sys.argv[2]
pattern = sys.argv[3]

# Open all of the files in the sample list as io.BufferedReader objects
tar = tarfile.open(tar_fn)
with open(sample_list) as samplefile:
    filenames = [pattern.format(s.strip()) for s in samplefile]
files = [tar.extractfile(fn) for fn in filenames]

n_files = len(files)

# Read the files in chunks, print them as you go, test if reached the end of file
while True:
    chunks = [f.read(CHUNKSIZE).decode("utf-8") for f in files]
    # f.read returns an empty bytestring on EOF, so break
    n_sites = len(chunks[0])
    if n_sites == 0:
        break
    else:
        for i_site in range(n_sites):
            site = ''.join(c[i_site] for c in chunks)
            if site[0]*n_files == site:
                print(site[0])
            else:
                print(site)

# Close all open files
for f in files:
    f.close()
tar.close()
