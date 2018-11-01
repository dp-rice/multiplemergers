mmc_genomics
------------

Investigating whether the Kingman coalescent can be rejected in favor of multiple mergers
    using genomic diversity data.

# Set-up
Clone repo `git clone [URL]`
## Conda
1. Install conda (https://conda.io/docs/index.html)
2. Create env: `conda env create -f config/environment.yml`
3. Activate environment: `source activate mmc`
## msprime-lambda (optional)
1. Clone repo: `git clone [URL]`
2. Install using `pip`: `cd [PATH]; pip install .`
## SLiM (optional)
Follow instructions at: https://messerlab.org/slim/
## fastNeutrino (optional)
https://sourceforge.net/projects/fastneutrino/

# Reproducing figures
Tarball with simulation output: `tar -xzf simulations.tgz`
files with processed DPGP3 data. `data/DPGP3/minor_allele_counts`
3 jupyter notebooks for making figures

# Reproducing analysis
run_snakemake.sh is for submitting to a SLURM cluster. (cluster.json)
2 snakefiles: one for running simulations, one for data analysis
Intructions on downloading raw dpgp3 data:
- http://www.johnpool.net/genomes.html
- `wget http://pooldata.genetics.wisc.edu/dpgp3_sequences.tar.bz2`
