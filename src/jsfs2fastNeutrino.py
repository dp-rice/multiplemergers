import sys

for line in sys.stdin:
    if line.startswith("#N_SAMPLES="):
        nsamples = line.split("=")[1]
        continue
    elif line.startswith("#"):
        continue
    sfs = line.split()

    print(nsamples, 1)
    print(100.0)
    for elem in sfs:
        print(elem)
    print(0.0)
