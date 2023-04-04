#!/bin/bash
#
# Two arguments: analysis name for reports and number of threads to use
#

#analysis="5tera"
analysis=$1

#threads=4
threads=$2

# map_jobs specifies how many alignment jobs to run concurrently; mapping jobs can be RAM exhaustive
map_jobs=1

date=$(date +"%Y%d%d_%H%M%S")

export TMPDIR=`pwd` # Tempdir for Snakemake 

CONDA_PATH=$(dirname $(dirname $(echo $(readlink -f `which conda`))))
source ${CONDA_PATH}/etc/profile.d/conda.sh

mkdir workflow/report

conda activate tera-snakemake

# Run 
snakemake \
    --use-singularity -c${threads} \
    --resources map_jobs=${map_jobs} \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.${analysis}.snakemake-stats.txt \
    -p

# add --rerun-incomplete if you run failed
# you might have to run snakemake --unlock snakemakedir as well

# Html report
snakemake \
    -c $threads \
    --report workflow/report/${date}.${analysis}snakemake-report.html 

# Text summary
snakemake \
    -c $threads \
    --detailed-summary > workflow/report/${date}.${analysis}.snakemake-summary.txt
