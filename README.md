

# Snakemake workflow: `TERA-Seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/opplatek/TERA-Seq_snakemake/workflows/Tests/badge.svg?branch=main)](https://github.com/opplatek/TERA-Seq_snakemake/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for `TERA-Seq data processing.`

This repository contains the main analysis steps used in the publication **TERA-Seq: True end-to-end sequencing of native RNA molecules for transcriptome characterization** (Ibrahim, F.\*, Oppelt, J.\*, Maragkakis, M.\* and Mourelatos, Z., 2021, Nucleic Acids Research, DOI: [10.1093/nar/gkab713](https://doi.org/10.1093/nar/gkab713), PMID: [34428294](https://pubmed.ncbi.nlm.nih.gov/34428294/). \* The authors wish it to be known that, in their opinion, the first three authors should be regarded as Joint First Authors. 

We kindly ask you to cite the publication above if you use any part of the analysis, code, or samples.

This is a **TERA-Seq Snakemake workflow**. It was designed to help you with TERA-Seq data analysis. As of now (04/08/2023), it is heavily based on the original workflow used in the publication (GitHub repository [TERA-Seq publication](https://github.com/mourelatos-lab/TERA-Seq_manuscript).  However, this is an **updated workflow** and certain steps don't have to be identical to the ones used in the original publication. 

## Usage

The usage of this workflow is described **in this README** and in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=opplatek%2FTERA-Seq_snakemake).

## Dependencies
This workflow is designed in Snakemake, and it is to be used with the TERA-Seq Singularity ([Docker container joppelt/teraseq:snakemake](https://hub.docker.com/repository/docker/joppelt/teraseq/general)) to ensure easy portability and reproducibility. 

There are only 3+1 main dependencies - **Git**, **Singularity**, **Snakemake**, and **Conda**. Conda is only used to install Snakemake and selected Python packages. You can install Snakemake independently on Conda. The extra Python packages are only used to create Snakemake reports and can also be omitted.  
### Git
If [Git](https://git-scm.com/) is not installed on your system (not very common), look for instructions on how to install it for your distro.

### Singularity
The workflow uses pre-installed software in a Singularity container. If you don't have Singularity installed on your system, you can get it from [here](https://docs.sylabs.io/guides/3.0/user-guide/index.html). We tested the workflow with Singularity `2.6.1-dist` on `openSUSE Leap 15.0` and `singularity-ce version 3.11.1-bionic` on `Ubuntu 18.04.6 LTS`.

We **do not** recommend installing Singularity using Conda as there might be permission issues.

#### Ubuntu 18.04 example
First, check whether you have all the dependencies from [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies) and [Go installed](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-go)

For `Ubuntu 18.04.6 LTS` and the tested version, do the following:
```bash
cd ~/tools/
wget https://github.com/sylabs/singularity/releases/download/v3.11.1/singularity-ce_3.11.1-bionic_amd64.deb
sudo apt install ./singularity-ce_3.11.1-bionic_amd64.deb && rm singularity-ce_3.11.1-bionic_amd64.deb
singularity version
``` 

### Snakemake

If you don't have Snakemake installed on your system, you can follow the instructions [here](https://snakemake.readthedocs.io/en/stable/index.html). The workflow was tested with Snakemake `7.24.0` and Python `3.7.1`. Although you can install Snakemake independently, we recommend using the provided YAML file to make a Conda environment. Please follow the instructions from the [next](#Conda) section.

### Conda
If you do not have Conda already installed, you can get it from [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). We tested the workflow on `conda 4.11.0` and `conda 4.12.0` installed through miniconda3.

The easiest way to install Snakemake is to use the provided YAML file and make the Conda environment. Before the Conda environment installation, we recommend installing [mamba](https://mamba.readthedocs.io/en/latest/installation.html) as it significantly speeds up the installation. However, it is not required. 

To install mamba with your **base** Conda environment activated:
```bash
conda install mamba -n base -c conda-forge
```

### Getting the workflow
Clone this GitHub repository to your chosen location. We recommend cloning the repository in your project directory to ensure proper reproducibility. 

Clone the latest commit with:
```bash
git clone https://github.com/mourelatos-lab/TERA-Seq_snakemake.git
```

You can install the Snakemake Conda environment as follows:
```bash
cd TERA-Seq_snakemake/
mamba env create -f environment.yaml -n teraseq-snakemake # If you don't have mamba, simply replace "mamba" at the beginning of the command with "conda"
```

If you don't want to use the provided YAML file, you can use the following:
```bash
mamba env create -n teraseq-snakemake \
    -c conda-forge -c bioconda -c anaconda -c r \
    python=3.7.1 snakemake-minimal=7.24.0 pandas pygments jinja2 networkx pygraphviz
```

To update the repository, you can execute the following:
```bash
cd TERA-Seq_snakemake/
git pull
```

#### Singularity container
The workflow **requires** Singularity container. Before making the environment, we recommend exporting temp and cache Singularity variables to a location with enough free disk space. To create one, please use the following:
```bash
tmpdir=$(pwd) # Your temporary and cache directory

export SINGULARITY_TMPDIR=$tmpdir
export SINGULARITY_CACHEDIR=$tmpdir

singularity pull docker://joppelt/teraseq:snakemake
```
This will create `teraseq-snakemake.simg` if you use Singularity 2 or  `teraseq_snakemake.sif` if you use Singularity 3.

## Workflow configuration
You must remember Snakemake is based on **relative** directory structure to the Snakefile or the workdir. This affects how we have to think about the config files and output directories.

### Config file
The main config file is in [`config/config.yaml`](config/config.yaml). The location of the main config file can be specified when launching the Snakemake, but we **strongly recommend** to keep it in your main project directory. 

You **have to** update the main config file with the correct paths to the sub-configs. Please note the paths can be either **absolute** or **relative** to Snakefile or workdir. We recommend copying the whole `config` directory to your workdir, keeping the sub-config files with relative paths, and only specifying the Singularity file directory name with an absolute path.

The main config files contain links to four sub-config files:
* **`samples.csv`**
A comma-delimited sample sheet with description and sample names to process. **Do not modify the header**. The sample sheet file has four columns:
    * `sample ID/name` - the name of the sample directory. The sample directory **must** have **`fastq`** subdirectory with **`reads.1.fastq.gz`** file inside. `reads.1.fastq.gz` are the basecalled reads created by the TERA-Seq protocol.
    * `assembly` - sample assembly. The implemented assemblies are: *hg38*, *mm10*, *sc3*.
    * `libtype` - TERA-Seq library type. The implemented library types are: *5tera*, *tera3*, *5tera3*.
    * `protocol` - RNA enrichment protocol. The implemented protocols are: *polya*, *total*.

For example: *`hsa.dRNASeq.HeLa.total.REL3.1,hg38,tera3,total`* is a *hsa.dRNASeq.HeLa.total.REL3.1* sample, the reference is *hg38*, it is *tera3* library, and *total* RNA was extracted.
* **`dirs.yaml `**
A YAML file specifying data, samples, and results directories. Please note these paths are **relative** to the Snakefile or workdir.
    * `datadir` - directory used to store references.
    * `samplesdir` - directory containing samples for the analysis. The samples inside this directory **must** be named **exactly** as in the `samples.csv` sample sheet. The individual samples **must** contain `fastq` directory with `reads.1.fastq.gz` FASTQ file.
    * `resdir` - directory name used to save the workflow results.
* **`sqldb`** 
Do you want to create an annotated SQL DB file? ["yes"|"no"]. **Must** include the double quotes. Default: "no".
* **`singularity_container`**
Absolute or relative path, including the name of the Singularity container file. This is most likely going to be `teraseq-snakemake.simg` for Singularity 2 or `teraseq_singularity.sif` for Singularity 3.
* **`resources.yaml`**
The maximum number of resources used throughout the workflow.
* **`ref_links.yaml`**
A YAML file specifying links to download the references. The workflow has been tested on [Ensembl](https://www.ensembl.org/index.html) references. 
This YAML file is structured as follows:
```yaml
organism reference name abbreviation:
    org: organism Latin name
    genome: link to download the reference genome
    gtf: link to download the reference genome annotation GTF file
    gtf_extend_fiveutr: number of nucleotides to extend the 5' UTR. New UTRs are created if they don't exist. Please note this has been tested only for the sc3 genome. Use "NA" if no nucleotides are to be added.
    gtf_extend_threeutr: number of nucleotides to extend the 3' UTR. New UTRs are created if they don't exist. Please note this has been tested only for the sc3 genome. Use "NA" if no nucleotides are to be added.
    gtrna: link to donwload [GtRNAdb](https://gtrnadb.ucsc.edu/) tRNA annotation archive.
    gtrna_bed: name of the decompressed BED file from the GtRNAdb tRNA annotation archive.
```
Use this file to add/update organisms/assemblies/references.

There is a special section for [Silva rRNA](https://www.arb-silva.de/) annotation:
```yaml
silva:
    lsu: Link to download Large (23S/28S) subunit Silva ribosomal RNAs FASTA file archive.
    ssu: Link to download Small (16S/18S) subunit Silva ribosomal RNAs FASTA file archive.
```
* **`adapters.yaml`**
YAML file with optimal adapter trimming parameters for individual TERA-Seq library types. Note: TERA-Seq has been tested on Nanopore Flow Cell ([R9.4.1 FLO-MIN106](https://store.nanoporetech.com/us/flow-cell-r9-4-1.html)) dRNA-Seq protocol ([SQK-RNA002](https://store.nanoporetech.com/us/direct-rna-sequencing-kit.html)). Optimal adapter removal settings will be different for other flow cells and protocols.

### Additional workflow resources
By default, we limit the maximum number of concurrent mapping jobs to four (you can change the default settings at [`config/resources.yaml`](config/resources.yaml) - look for `map_jobs`). There are two reasons - a) classical HDD might have trouble running more than ~4 I/O demanding jobs simultaneously; b) `sam-count-secondary` can consume a lot of RAM per sample (depending on the sequencing depth). If you run the analysis on SSD and/or have a lot of RAM, you can increase the maximum value or add `--resources map_jobs=maxnumberofmapjobs` to the Snakemake command (replace `maxnumberofmapjobs` with a value of the maximum number of concurrent mapping jobs).

## Running the workflow
Make sure you activate your Conda environment if this is how you installed Snakemake:
```bash
conda activate teraseq-snakemake
```
The workflow is designed to run in any work directory. 

Let's assume you cloned the main TERA-Seq_snakemake repository to `/home/user/tools/TERA-Seq_snakemake`, and pulled the Singularity container to the same directory. Your main project directory is `/home/user/projects/TERASeq`. The main project directory contains individual sample directories in `data/samples`. Each of the sample's directories has `fastq/reads.1.fastq.gz` file. 

1. Copy the `config` directory from the main TERA-Seq_snakemake repository copy:
```bash
cp -r /home/user/tools/TERA-Seq_snakemake/config /home/user/projects/teraseq/
```
2. Modify the main config to reflect the location of the Singularity container (**absolute path**) set if you also want to make the SQLite db file. As stated previously, the main config file can be found at `config/config.yaml` (See [Config file section](#config-file) for more details). After the modifications, the config might look like this:
```yaml
samples: "config/samples.csv" # Sample sheet

dirs: "config/dirs.yaml" # Name of directories for samples, references, and results

sqldb: "yes" # [yes|no] Make SQL db file

singularity_container: "/home/user/tools/TERA-Seq_snakemake/teraseq_snakemake.sif" # If we are using Singularity 3
#singularity_container: "/home/user/tools/TERA-Seq_snakemake/teraseq_snakemake.simg" # If we are using Singularity 2

ref_links: "config/ref_links.yaml" # Links to download reference

adapters: "config/adapters.yaml" # Adapter trimming settings
```
3. Modify the samples, references, and results output directories if necessary at `config/dirs.yaml`. See [Config file section](#config-file) for more details.
4. Add the samples to analyze to `config/samples.csv`. See [Config file section](#config-file) for more details.
5. Run a dry-run to test we set everything correctly:
```bash
conda activate teraseq-snakemake

threads=16 
concurrent_mappings=4 

snakemake \
    -c $threads --use-singularity \
    --resources map_jobs=$concurrent_mappings \
    --snakefile /home/user/tools/TERA-Seq_snakemake/workflow/Snakefile \
    --directory /home/user/projects/TERASeq \
    --configfile config/config.yaml \
    -pn
```
6. Run the pipeline with workflow stats and reports (Note: the reports need the additional Python packages mentioned in the [Getting the worflow section](getting-the-workflow):
```bash
date=$(date +"%Y%d%d_%H%M%S")
mkdir report

# Run the workflow with stats
snakemake \
    -c $threads --use-singularity \
    --resources map_jobs=$concurrent_mappings \
    --snakefile /home/user/tools/TERA-Seq_snakemake/workflow/Snakefile \
    --directory /home/user/projects/TERASeq \
    --configfile config/config.yaml \
    --stats report/${date}.myFirstRun.teraseq-snakemake-stats.txt
    -p

# Make reports and summaries
snakemake \
    -c $threads \
    --report report/${date}.myFirstRun.teraseq-snakemake-report.html 

snakemake \
    -c $threads \
    --detailed-summary > report/${date}.myFirstRun.teraseq-snakemake-summary.txt
```

### Helper run script
If you don't want to remember how to run Snakemake, you can use the *helper* [`run.sh`](run.sh) script. Run:
```bash
# Use the same settings as in the example
./run.sh \
    /home/user/tools/TERA-Seq_snakemake/workflow/Snakefile \
    /home/user/projects/TERASeq \
    config/config.yaml \
    "myFirstRun" \
    16 \
    4

# Use the default settings only with mandatory arguments
./run.sh /home/user/tools/TERA-Seq_snakemake/workflow/Snakefile /home/user/projects/TERASeq 
```

## Results
### 5TERA sample
```bash
5tera-sample-ID
|-- fastq
|   |-- reads.1.fastq.gz # Original FASTQ file
|   |-- reads.1.sanitize.fastq.gz # Read name sanitized FASTQ file
|   |-- reads.1.sanitize.adapt_trim.5tera.fastq.gz # 5TERA adapter-removed FASTQ file - source
|   |-- reads.1.sanitize.w_adapt.5tera.fastq.gz # FASTQ file with reads w/ 5TERA adapter found (adapter removed)
|   |-- reads.1.sanitize.w_adapt.5tera.names.txt # List of reads w/ 5TERA adapter found
|   |-- reads.1.sanitize.wo_adapt.5tera.fastq.gz # FASTQ file with reads w/o 5TERA adapter found
|   |-- reads.1.sanitize.wo_adapt.5tera.names.txt # List of reads w/o 5TERA adapter found
|   |-- reads.1.sanitize.adapt_trim.fastq.gz -> reads.1.sanitize.adapt_trim.5tera.fastq.gz # 5TERA adapter-removed FASTQ file - link
|   `-- reads.1.sanitize.noribo.fastq.gz # FASTQ w/o rRNA mapped reads
|-- align
|   |-- reads.1.sanitize.toGenome.sorted.bam # Genome mapping BAM
|   |-- reads.1.sanitize.toGenome.sorted.bam.bai # Genome mapping BAM index
|   |-- reads.1.sanitize.toRibosomal.sorted.bam # rRNA-only mapping BAM
|   |-- reads.1.sanitize.noribo.toTranscriptome.sorted.bam -> reads.1.sanitize.noribo.toTranscriptome.trans.sorted.bam # Transcriptome mapping BAM - link
|   `-- reads.1.sanitize.noribo.toTranscriptome.trans.sorted.bam # Transcriptome mapping BAM - source
|-- db # Only created if sqldb=="yes"
|   |-- *.done # Various SQL DB annotation checkpoint files
|   `-- sqlite.db # SQL DB annotated file with both genome and transcriptome mappings
`- log
    |-- cutadapt.5tera.log # 5TERA adapter Cutadapt log file
    |-- cutadapt.len.tsv # Length of removed adapters
    |-- reads.1.sanitize.adapt_trim.read-len.tsv.gz # Adapter-removed read lengths
    |-- cutadapt.len.pdf # Visualization of removed adapter length
    |-- reads.1.sanitize.toGenome.align-len.full.tsv.gz # Aligned portion of read lengths (excludes soft clipping and similar) - genome mapping/alignment stats
    |-- reads.1.sanitize.toGenome.mapped-len.full.tsv.gz # Mapped read lengths (includes soft clipping and similar) - genome mapping/alignment stats
    |-- reads.1.sanitize.toGenome.len.mapped-vs-aligned.pdf # Visualization of mapped read length and the actual aligned portion of the read - genome mapping/alignment stats
    |-- reads.1.sanitize.toGenome.len.trim-vs-aligned.pdf # Visualization of adapter-removed read length and the actual aligned portion of the read - genome mapping/alignment stats 
    |-- reads.1.sanitize.noribo.toTranscriptome.align-len.full.tsv.gz # Aligned portion of read lengths (excludes soft clipping and similar) - transcriptome mapping/alignment stats 
    |-- reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.tsv.gz # Mapped read lengths (includes soft clipping and similar) - transcriptome mapping/alignment stats
    |-- reads.1.sanitize.noribo.toTranscriptome.len.mapped-vs-aligned.pdf # Visualization of mapped read length and the actual aligned portion of the read - transcriptome mapping/alignment stats
    |-- reads.1.sanitize.noribo.toTranscriptome.len.trim-vs-aligned.pdf # Visualization of adapter-removed read length and the actual aligned portion of the read - transcriptome mapping/alignment stats 
    |-- reads.1.sanitize.toGenome.mapping-stats.tsv # Basic (samtools stat) mapping statistics - genome stats
    |-- sqldb.mapping-stats.transcriptome.txt # Basic (sqldb-based) mapping statistics - transcriptome stats; only created if sqldb=="yes"
    |-- sqldb.mapping-stats.genome.txt # Basic (sqldb-based) mapping statistics - genome stats; only created if sqldb=="yes"
    `-- reads.1.sanitize.noribo.toTranscriptome.mapping-stats.tsv # Basic (samtools stat) mapping statistics - transcriptome stats
```

### TERA3 sample - differences from 5TERA sample only
```bash
|-- fastq
|   |-- reads.1.sanitize.adapt_trim.fastq.gz -> reads.1.sanitize.adapt_trim.tera3.fastq.gz # TERA3 adapter-removed FASTQ file - link
|   |-- reads.1.sanitize.adapt_trim.tera3.fastq.gz # TERA3 adapter-removed FASTQ file - source
|   |-- reads.1.sanitize.w_adapt.tera3.fastq.gz # FASTQ file with reads w/ TERA3 adapter found (adapter removed)
|    `-- ...
|   |-- reads.1.sanitize.w_adapt.tera3.names.txt # List of reads w/ TERA3 adapter found
|   |-- reads.1.sanitize.wo_adapt.tera3.fastq.gz # FASTQ file with reads w/o TERA3 adapter found
|   |-- reads.1.sanitize.wo_adapt.tera3.names.txt # List of reads w/o TERA3 adapter found
|    `-- ...
|-- align # Identical to the 5TERA example
|    `-- ...
|-- db # Identical to the 5TERA example
|    `-- ...
`-- log
    |-- cutadapt.tera3.full-length.log # TERA3 adapter Cutadapt log file from the whole reads after subsetting only for TERA3-found reads
    |-- cutadapt.tera3.log # TERA3 adapter Cutadapt log file only from the last 200 nt
    `-- ...
```

### 5TERA3 sample - differences from 5TERA and TERA3 samples only
5TERA3 results are a combination of 5TERA and TERA3 results. There are no unique files.

## TODO list
* Add automatic detection of singularity container in case the path is not set in the config file
* Check for 5utr and 3utr extensions in the references config list as automatically set to NA if not set in the config
* Expand initial QC tests
* Implement analyses from the TERA-Seq manuscript

> Written with [StackEdit](https://stackedit.io/).
