#!python3

import gzip as gz
import os


# ----- Definitions for the project ----- #

# CHROMS = [i for i in range(22,0, -1)]

#REF_GENOME = "/home/abiddanda/novembre_lab/data/external_public/reference_genomes/hs37d5.fa"

#PBWT = 'bin/pbwt/pbwt'

#DATA_PATH = "/home/abiddanda/novembre_lab/share/botai/data/"

#GEODIST_FILE = "/home/abiddanda/novembre_lab/data/external_public/geodist/newll"

#GEODIST_PROJECT = "/home/abiddanda/novembre_lab/abiddanda/hackathonDec2016"

#INDIVS = ['BOT2016', 'Yamnaya', 'Yana', 'Kolyma_River']

# ---- Useful Functions ---- #

#base = lambda x:os.path.splitext(x)[0]

# ---- Including other snakefiles ---- #
include: 'snakefiles/dpgp3.snake'
include: 'snakefiles/simulations.snake'

# ---- Base Rule ----------- #

rule targets:
    input:
        "sandbox/test.txt"

rule test:
    output:
        "sandbox/test.txt"
    shell:
        "echo Test > {output}"
