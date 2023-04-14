#!/bin/bash
#
# Five arguments: Snakefile path (required), Analysis/workdir path (required), analysis name for the reports (optional, default: "run"), and number of threads to use (optional, default: 1), and number of concurrent mapping jobs (optional, default: 4)
#
################################
usage="Usage: $(basename "$0") analysisName numberOfThreads (default: 1)"

if [ "$#" -eq 0 ] | [ "$#" -lt 2 ]; then
    echo $usage
    exit 1
fi

snakefile=$1

workdir=$2

# Analysis name (for reports)
analysis="run"
if [ "$#" -eq 3 ]; then
    analysis=$3
fi

# Number of threads to use
threads=1
if [ "$#" -eq 4 ]; then
    threads=$4
fi

# The maximum number of concurrent mapping jobs
map_jobs=4
if [ "$#" -eq 5 ]; then
    map_jobs=$5
fi

### Additional variables
# map_jobs specifies how many alignment jobs to run concurrently; mapping jobs can be RAM exhaustive
export TMPDIR=`pwd` # Tempdir for Snakemake
date=$(date +"%Y%d%d_%H%M%S")
reportdir="workflow/report"

################################

echo "Running >>> ${analysis} <<< analysis with ${snakefile} Snakefile in ${workdir} workdir with ${threads} threads, and with the maximum of ${map_jobs} mapping jobs."
 
################################
# In case Conda doesn't load automatically
CONDA_PATH=$(dirname $(dirname $(echo $(readlink -f `which conda`))))
source ${CONDA_PATH}/etc/profile.d/conda.sh

conda activate teraseq-snakemake

mkdir -p $reportdir

# Run 
snakemake \
    --use-singularity -c ${threads} \
    --resources map_jobs=${map_jobs} \
    --snakefile ${snakefile} \
    --directory ${workdir} \
    --configfile config/config.yaml \
    --stats ${reportdir}/${date}.${analysis}.teraseq-snakemake-stats.txt \
    -p

# Add --rerun-incomplete if your run failed and didn't finish successfully
# You might need to run 'snakemake --unlock snakemakedir' as well

# Html report
snakemake \
    -c ${threads} \
    --report ${reportdir}/${date}.${analysis}.teraseq-snakemake-report.html 

# Text summary
snakemake \
    -c ${threads} \
    --detailed-summary > ${reportdir}/${date}.${analysis}.teraseq-snakemake-summary.txt
