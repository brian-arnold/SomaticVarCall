#!/bin/bash

source /n/fs/ragr-research/users/bjarnold/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
snakemake --cores 50 --use-conda
#snakemake -n -p 


