#!/bin/bash
#
# Two arguments: analysis name for reports and number of threads to use
#
# Note: By default, we limit number of alignment jobs running concurrently to 1. If you have enough RAM 
#           (the most limiting factor) you can increase map_jobs variable
#
################################
usage="Usage: $(basename "$0") analysisName numberOfThreads (defaul: 1)"

if [ "$#" -eq 0 ] | [ "$#" -lt 2 ]; then
    echo $usage
    exit 1
fi

# Analysis name (for reports)
#analysis="5tera"
analysis=$1

# Number of threads to use
threads=1
if [ "$#" -eq 2 ]; then
    threads=$2
fi

### Additional variables
# map_jobs specifies how many alignment jobs to run concurrently; mapping jobs can be RAM exhaustive
map_jobs=1
export TMPDIR=`pwd` # Tempdir for Snakemake
date=$(date +"%Y%d%d_%H%M%S")
reportdir="workflow/report"

echo "Running >>> ${analysis} <<< analysis with ${threads} threads with maximum of ${map_jobs} mapping jobs."
 
################################
CONDA_PATH=$(dirname $(dirname $(echo $(readlink -f `which conda`))))
source ${CONDA_PATH}/etc/profile.d/conda.sh

conda activate teraseq-snakemake
################################

mkdir $reportdir

# Run 
snakemake \
    --use-singularity -c${threads} \
    --resources map_jobs=${map_jobs} \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.${analysis}.snakemake-stats.txt \
    -p

# Add --rerun-incomplete if your run failed and didn't finish successfully
# You might need to run 'snakemake --unlock snakemakedir' as well

# Html report
snakemake \
    -c $threads \
    --report workflow/report/${date}.${analysis}snakemake-report.html 

# Text summary
snakemake \
    -c $threads \
    --detailed-summary > workflow/report/${date}.${analysis}.snakemake-summary.txt
