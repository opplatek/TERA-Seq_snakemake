#!/bin/bash
#
# Six arguments: 
#   Snakefile path (default: workflow/Snakefile)
#   Analysis/workdir path (default: path to this script)
#   Main config file path (default: config/config.yaml)
#   Analysis name for the reports (default: "run")
#   Number of threads to use (default: 1)N
#   Number of concurrent mapping jobs (default: 4)
#
################################
usage="Usage: $(basename "$0") snakefilePath (default: workflow/Snakefile) workDir (default: path to this script) configPath (default: config/config.yaml) analysisName (default: \"run\") numberOfThreads (default: 1) numerOfMapJobs (default: 4)"

if [ "$#" -eq 0 ]; then
    echo $usage
    exit 1
fi

# Snakefile path
snakefile="workflow/Snakefile"
if [ "$#" -ge 1 ]; then
    snakefile=$1
fi

workdir=$(dirname $(realpath -s $0))
if [ "$#" -ge 2 ]; then
    workdir=$2
fi

# Main config file
config="config/config.yaml"
if [ "$#" -eq 3 ]; then
    config=$3
fi

# Analysis name (for reports)
analysis="run"
if [ "$#" -eq 4 ]; then
    analysis=$4
fi

# Number of threads to use
threads=1
if [ "$#" -eq 5 ]; then
    threads=$5
fi

# The maximum number of concurrent mapping jobs
map_jobs=4
if [ "$#" -eq 6 ]; then
    map_jobs=$6
fi

### Additional variables
# map_jobs specifies how many alignment jobs to run concurrently; mapping jobs can be RAM exhaustive
export TMPDIR=`pwd` # Tempdir for Snakemake
date=$(date +"%Y%d%d_%H%M%S")
reportdir="workflow/report"

################################

echo "Running >>> ${analysis} <<< analysis with ${snakefile} Snakefile in ${workdir} workdir, with main config ${config}, ${threads} threads, and a maximum of ${map_jobs} mapping jobs."
 
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
    --configfile ${config} \
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
