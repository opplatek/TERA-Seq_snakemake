# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `<description>`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.# Singularity and Snakemake

## Notes
* Nice summary [here](https://siscourses.ethz.ch/container_pipeline_tutorial/2_container_pipeline_tutorial.pdf)
* Another decent summary [here](https://wfbroderick.com/2022-Aug-01.html) and [here](https://snakemake-on-nesi.sschmeier.com/singularity.html)
* General simple Docker build info [here](https://devopscube.com/build-docker-image/)

## *Inside* vs *outside* 
* Running Singularity has two "main" options
    1) Bind the data into the Singularity container and run the whole analysis inside - less flexibility but more reproducibility (cannot use host-installed software)
    2) Specify to run Singularity in each Snakemake rule and collect the results - more flexibility but less reproducibility (needs Python, Snakemake, ..., installed on the host)

## (?) Our solution
* Install Python3, Conda (Miniconda), and Snakemake (inside Conda) on the host and only run the individual steps in the Singularity container
    * We'll hope most of the hosts have Python3 installed and Miniconda doesn't need root access
* Or pre-make Singularity container like [here](https://snakemake-on-nesi.sschmeier.com/singularity.html) and run the whole workflow with `snakemake --use-singularity`
    
## Install Conda environment
```
cd /home/jan/projects/TERA-Seq_snakemake/
conda install mamba -n base -c conda-forge
mamba env create -f environment.yaml -n tera-snakemake # If you don't have mamba, replace it with conda
# conda env create -f environment.yaml -n tera-snakemake
```

### Install Singularity
**Important:** Installing Singularity in Conda doesn't work well. The system doesn't see it as a root installation and converts all Singularity containers (.sif) into a sandbox for **every** Snakemake output.

First, check whether you have all the dependencies from [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies) and have [Go installed](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-go)

For Ubuntu 20.04 and the tested version, do:
```
cd ~/tools/
wget https://github.com/sylabs/singularity/releases/download/v3.11.1/singularity-ce_3.11.1-bionic_amd64.deb
sudo apt install ./singularity-ce_3.11.1-bionic_amd64.deb && rm singularity-ce_3.11.1-bionic_amd64.deb
singularity version
```

### Make singularity image
singularity pull docker://joppelt/teraseq-snakemake

### Run the pipeline
```
conda activate tera-snakemake
snakemake --use-singularity -c1 -p
```
  
## TODO
* How to run different Conda environments within the Singularity?
    * I like to pre-generate all Conda envs Snakemake workflow needs - can this be done in Singularity? I guess you would have to activate Conda on host, and then have pre-made Conda envs in the container and activate them in (?) Snakemake rule wrappers
* Cut Dockerfile to the bare minimum (if it works)
* If the *old* Dockerfile doesn't work, make a new one

### Make Singularity from Docker
```
#docker build -t local/my_container:latest .
sudo docker build -t teraseq:snakemake .
#sudo singularity build my_container.sif docker-daemon://local/my_container:latest
mkdir /home/jan/tmp # Make temporary directory somewhere fast
export SINGULARITY_CACHEDIR=/home/jan/tmp # Export Singularity cache variable somewhere fast otherwise it will jam your /tmp
singularity build --tmpdir '/home/jan/tmp' teraseq-snakemake.sif docker-daemon://local/teraseq:snakemake # I haven't figured out how to use variable in Singularity --tmpdir so you have to manually write the path
'''

* run Singularity in an interactive mode
`singularity shell teraseq-snakemake.sif`

## Notes
* `docker run -ti my_image /bin/bash`
* `docker build -t local/teraseq:perl - < Dockerfile-perlOnly`
* nice workflow [here](https://github.com/zavolanlab/zarp)
