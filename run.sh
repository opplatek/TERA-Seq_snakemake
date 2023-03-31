#!/bin/bash
#
# Two arguments: analysis name for reports and number of threads to use
#

#analysis="5tera"
analysis=$1

#threads=4
threads=$2

date=$(date +"%Y%d%d_%H%M%S")

export TMPDIR=`pwd` # Tempdir for Snakemake 

CONDA_PATH=$(dirname $(dirname $(echo $(readlink -f `which conda`))))
source ${CONDA_PATH}/etc/profile.d/conda.sh

mkdir workflow/report

conda activate tera-snakemake

# Run 
snakemake \
    --use-singularity -c${threads} \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.${analysis}.snakemake-stats.txt \
    -p

# Html report
snakemake \
    -c $threads \
    --report workflow/report/${date}.${analysis}snakemake-report.html 

# Text summary
snakemake \
    -c $threads \
    --detailed-summary > workflow/report/${date}.${analysis}.snakemake-summary.txt
