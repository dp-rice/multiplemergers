#!python3

import gzip as gz
import os


# ----- Definitions for the project ----- #

#CHROMS = [i for i in range(22,0, -1)]

#REF_GENOME = "/home/abiddanda/novembre_lab/data/external_public/reference_genomes/hs37d5.fa"

#PBWT = 'bin/pbwt/pbwt'

#DATA_PATH = "/home/abiddanda/novembre_lab/share/botai/data/"

#GEODIST_FILE = "/home/abiddanda/novembre_lab/data/external_public/geodist/newll"

#GEODIST_PROJECT = "/home/abiddanda/novembre_lab/abiddanda/hackathonDec2016"

#INDIVS = ['BOT2016', 'Yamnaya', 'Yana', 'Kolyma_River']

# ---- Useful Functions ---- #

#base = lambda x:os.path.splitext(x)[0]

# ---- Including other snakefiles ---- #
#include: 'snakefiles/sfs.snake'
#include: 'snakefiles/geodist.snake'
#include: 'snakefiles/a_share.snake'


# ---- Base Rule ----------- #

rule targets:
    input:
        "sandbox/test.txt"

rule test:
    output:
        "sandbox/test.txt"
    shell:
        "echo Test > {output}"
    

#rule all:
#	input:		
#            expand('data/1kg_phase3_sfs/ALL.chr{CHROM}.{pop}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.frq.count', pop=POPS, CHROM=[22]),
#            expand('plots/geodist/{name}_maf_{maf}_total_1kg_geodist3.png', name=INDIVS, maf=[0.01, 0.05]),
#            expand('data/logl_ashare/total/{maf}.logl.hwe_table.txt', maf=[0.1, 0.05, 0.02])
