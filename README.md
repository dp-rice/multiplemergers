# Distinguishing multiple-merger from Kingman coalescence using two-site frequency spectra

### Daniel P. Rice, John Novembre, and Michael M. Desai

This repository accompanies [Rice, Novembre, and Desai (2018)](https://www.biorxiv.org/content/early/2018/11/03/461517). It includes:

1. `jupyter` notebooks, along with processed data and simulation output, for reproducing the figures in the manuscript.
2. `snakemake` pipelines for reproducing simulations and data analysis from scratch.
3. `python` source code for computing two-site frequency spectra (2-SFS) and frequency pointwise mutual information (fPMI) in Kingman coalescent models of population growth and beta coalescent models.

# Getting started
First, clone this repository using:

```git clone https://github.com/dp-rice/multiplemergers.git```

## Setting up the computing environment

The repository contains a `conda` environment file for reproducing the computing environment where the project was originally carried out. Conda will make sure that your `python` version is the same one we used and that you have all of the necessary packages, including `jupyter` and `snakemake`. To use this environment:

1. Install `conda`, following the instructions [here](https://conda.io/docs/index.html).

2. From inside the project directory, create a new `conda` environment using:

 `conda env create -f config/environment.yml`

3. Activate the `conda` environment:

`source activate mmc`

See the `conda` documentation for more information about using environments.

## Coalescent simulations using `msprime-lambda` (optional)

The coalescent simulations were performed using a modified version of `msprime`, which allows for multiple mergers (lambda) coalescence.
Most of the under-the-hood work was done by [Joe Zhu](https://github.com/shajoezhu), but we made some modifications to allow access from the python API.
We would intend to contribute this feature to the main `msprime` codebase.

For now, if you would like to reproduce our multiple mergers coalescent simulations, you will need to install our version of `msprime`:

1. Clone the repository wherever you would like on your machine:

 `git clone https://github.com/dp-rice/msprime-lambda.git`

2. With the `conda` environment activated, install `msprime` locally using `pip`:

`cd [BASE/DIRECTORY/OF/MSPRIME-LAMBDA]; pip install .`

You will now be able to use `import msprime` within `python` code to use the msprime python API whenever the `mmc` conda environment is activated.

## Forward-time simulations using `SLiM` (optional)

In order to reproduce the selective sweep simulations, you will need to install `SLiM`.
Follow the instructions at: https://messerlab.org/slim/

## Fitting demographic models with `fastNeutrino` (optional)
In order to reproduce our coalescent model fits to *Drosophila melanogaster* diversity data, you will need to install `fastNeutrino`.
Follow the instructions at: https://sourceforge.net/projects/fastneutrino/

# Reproducing figures
If you would like to reproduce the figures in the manuscript without re-running the simulations or data pre-processing, we have provided the output of these steps.

Our simulation output is in a `tar` archive. To unpack the archive, use the command: `tar -xzf simulations.tgz`.

The processed *Drosophila* data is in the directory `data/DPGP3/minor_allele_counts`. Each file contains the number of genotyped samples (out of 100) and minor allele count at every site in one chromosome arm.

The commands for reproducing our figures are in three `jupyter` notebooks:

1. `notebooks/figures.ipynb` contains computes two-site frequency spectra from theoretical calculations and simulation output and generates most of the figures in the paper.

2. `notebooks/figure5.ipynb` produces Figure 5 (hiloPMI vs. genetic distance) from simulation output.

3. `notebooks/DPGP3_figures.ipynb` computes the 2-SFS from *Drosophila* data and compares it to coalescent simulations of population growth.

# Reproducing simulations and data processing
If you would like to reproduce the analysis from scratch or do something similar, we have provided two `snakemake` pipelines:

1. `snakefiles/simulations.snake` contains the pipeline for running all of the simulations.

2. `snakefiles/dpgp3.snake` contains the pipeline for processing the *Drosophila* data. You will need to download the raw data from the [Drosophila Genome Nexus](http://www.johnpool.net/genomes.html). You can also find a description of the data there. To download the data, use:

`wget http://pooldata.genetics.wisc.edu/dpgp3_sequences.tar.bz2`

To run snakemake locally, use commands like `snakemake msprime_all`, which will run all the `msprime` simulations. See the snakefiles for all rules. There is also a bash script `run_snakemake.sh`, which will submit jobs to a parallel computing cluster. It is currently configured for `SLURM` job scheduling.

# Calculating the 2-SFS and frequency PMI
The `src` directory contains `python` code for computing the various statistics described in the manuscript.

- `src/coalescentmoments.py` contains functions for calculating the 2-SFS for the beta coalescent.
- `src/zivkovic.py` contains functions for calculating the 2-SFS for the Kingman coalescent with population growth.
- `src/simulate_joint_sfs.py` is a command-line interface for running msprime simulations and calculating the 2-SFS (sometimes called "joint SFS in the code").

For examples of how to use these functions, see the notebooks and `snakemake` files.
