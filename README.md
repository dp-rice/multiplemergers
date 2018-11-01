mmc_genomics
------------

Investigating whether the Kingman coalescent can be rejected in favor of multiple mergers
    using genomic diversity data.

dependencies:
- conda
- conda environment
- msprime-lambda: git clone [REPO]; cd [DIRECTORY]; pip install .
- SLiM
- fastNeutrino

run_snakemake.sh is for submitting to a SLURM cluster. (cluster.json)
2 snakefiles: one for running simulations, one for data analysis
3 jupyter notebooks for making figures
Tarball with simulation output
Intructions on downloading raw dpgp3 data.
files with processed DPGP3 data.
