import sys

script, n_old, n_new = sys.argv
n_old = int(n_old)
n_new = int(n_new)

infile = sys.stdin
outfile = sys.stdout

# If the two sample sizes are equal, copy the file
if n_old == n_new:
    for line in infile:
        outfile.write(line)
    exit()


line = infile.readline() 
# Re-print header
while not line.startswith('#OUT'):
    outfile.write(line)
    line = infile.readline()

# Print the correct sample size
sline = line.split()
sline[-1] = str(n_new)
outfile.write(' '.join(sline) + '\n')
# This is the line that says 'Mutations:'
outfile.write(infile.readline())

# Get mutations
mutations = {}
while True:
    line = infile.readline()
    if line.startswith('Genomes:'):
        break
    else:
        sline = line.split()
        # Change the whole-sample prevalence to a counter starting at zero
        sline[-1] = 0
        mutations[sline[0]] = sline

# Get first n_new samples
genomes = [infile.readline().split() for i in range(n_new)]

# Count mutations
for gen in genomes:
    # The mutation list starts at position 2
    for mut_id in gen[2:]:
        # Increment the counter
        mutations[mut_id][-1] += 1

# Write mutations
for mut in mutations:
    sline = mutations[mut]
    counts = sline[-1]
    if counts > 0 and counts < n_new:
        # Replace whole-sample prevalence with subsample prevalence
        sline[-1] = str(counts)
        outfile.write(' '.join(sline) + '\n')

# Write subsampled genomes
outfile.write('Genomes:\n')
for gen in genomes:
    outfile.write(' '.join(gen) + '\n')
