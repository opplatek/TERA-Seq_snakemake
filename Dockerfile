# Download base image ubuntu 16.04
FROM ubuntu:16.04

# LABEL about the custom image
LABEL maintainer="jan.oppelt@pennmedicine.upenn.edu"
LABEL version="0.1"
LABEL description="This is custom Docker Image for \
analysis of TERA-Seq publication (DOI: https://doi.org/10.1093/nar/gkab713)."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Set default shell
SHELL ["/bin/bash", "-c"]

### System-wide requirements; cpanminus is not required if Perl uses virtual environment method; g++ and zlib1g-dev are required only for Nanopolish
RUN apt-get update \
    && apt-get install -y git gcc make wget g++ zlib1g-dev cpanminus curl locales language-pack-en \
    && locale-gen en_US en_US.UTF-8 && dpkg-reconfigure locales \
    && rm -rf /var/lib/apt/lists/*

### Set locale settings
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV LC_TIME en_US.UTF-8

### Main GitHub repo
WORKDIR /usr/local
RUN git clone https://github.com/mourelatos-lab/TERA-Seq_manuscript.git

### Install Miniconda3
ENV PATH "/usr/local/miniconda3/bin:${PATH}"
ARG PATH="/usr/local/miniconda3/bin:${PATH}"
RUN echo ${PATH}

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O Miniconda3.sh \
    && mkdir /usr/local/.conda \
    && bash Miniconda3.sh -bfp /usr/local/miniconda3 \
    && rm -f Miniconda3.sh
#RUN conda --version

## Install Mamba for faster installation
#RUN conda install -c conda-forge mamba

# Get Conda yml and install environment
#RUN mamba env create -f /usr/local/TERA-Seq_manuscript/teraseq-env.yml
RUN conda env create -f /usr/local/TERA-Seq_manuscript/teraseq-env.yml

### Fix error when using preinstalled Conda envs activation inside a Singularity container (https://gitlab.univ-nantes.fr/bird_pipeline_registry/srp-pipeline/-/tree/b5c1ac5e4f2449605701040484b0905028a74767#Troubleshooting or https://github.com/conda/conda/issues/8186)
#       /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-binutils_linux-64.sh: line 67: HOST: unbound variable
# Start Signularity in an interactive (singularity shell container.sf) mode and check what's the problem (sed -n '67p') and follow the same pattern:
#     "You can fix it by adapting the previous trick on the problematic lines (e.g. line 67 eval oldval="\$${from}$thing" becomes eval oldval="\${${from}$thing:-}" and so on)"
# You might get a similar error with virtual Perl environment. In that case, https://github.com/conda/conda/issues/8186#issuecomment-532874667 solution was enough - in the 'script', add set +eu, load the environment, and set set -eu
RUN sed -i 's/eval\ oldval="\\$\${from}\$thing"/eval oldval="\\${${from}$thing:-}"/' /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-binutils_linux-64.sh \
    && sed -i 's/\${SYS_SYSROOT}/${SYS_SYSROOT:-}/' /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-gcc_linux-64.sh \
    && sed -i 's/eval\ oldval="\\$\${from}\$thing"/eval oldval="\\${${from}$thing:-}"/' /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-gcc_linux-64.sh \
    && sed -i 's/eval\ oldval="\\$\${from}\$thing"/eval oldval="\\${${from}$thing:-}"/' /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-gfortran_linux-64.sh \
    && sed -i 's/eval\ oldval="\\$\${from}\$thing"/eval oldval="\\${${from}$thing:-}"/' /usr/local/miniconda3/envs/teraseq/etc/conda/activate.d/activate-gxx_linux-64.sh

# Increase default FastQC RAM
RUN sed -i 's/-Xmx250m/-Xmx5g/g' /usr/local/miniconda3/envs/teraseq/opt/fastqc-*/fastqc

#ENV PATH="${PATH}:/usr/local/miniconda3/envs/teraseq/bin"

RUN ln -s /usr/local/miniconda3/envs/teraseq/bin/R /bin/R
#   \ && ln -s /usr/local/miniconda3/envs/teraseq/bin/curl /bin/curl

### Save default Conda path
RUN sed -i '/CONDA_PREFIX/d' /usr/local/TERA-Seq_manuscript/PARAMS.sh \
    && echo -e "CONDA_PREFIX=\"/usr/local/miniconda3\"" >> /usr/local/TERA-Seq_manuscript/PARAMS.sh

# Set permissions the same for all users to avoid Singularity "issue" with running as root (Docker) vs running as user who executed the command (Singularity)
RUN chmod -R a=u /usr/local/miniconda3

### Perl
WORKDIR /usr/local/TERA-Seq_manuscript/tools

## Virtual environment install (option 2)
RUN git clone https://github.com/jizhang/perl-virtualenv.git \
    && cd perl-virtualenv/ \
    && git reset f931774 --hard \
    && chmod u+x virtualenv.pl \
    && ./virtualenv.pl teraseq \
    && . teraseq/bin/activate \
    && curl -L https://cpanmin.us/ -o teraseq/bin/cpanm \
    && chmod +x teraseq/bin/cpanm

## TODO: Add versions to all Perl modules installations
RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm inc::Module::Install \
    && cpanm autodie \
    && cpanm DBI \
    && cpanm Devel::Size \
    && cpanm Getopt::Long::Descriptive \
    && cpanm IO::File \
    && cpanm IO::Interactive \
    && cpanm IO::Uncompress::Gunzip \
    && cpanm Params::Validate \
    && cpanm Params::Util \
    && cpanm Sub::Install \
    && cpanm Modern::Perl \
    && cpanm --force MooseX::App::Simple \
    && cpanm --force MooseX::App::Command \
    && cpanm --force MooseX::Getopt::Meta::Attribute::Trait::NoGetopt

RUN git clone --recursive https://github.com/genoo/GenOO.git perl-virtualenv/teraseq/lib/perl5/GenOO_git \
    && cd perl-virtualenv/teraseq/lib/perl5/GenOO_git/ \
    && git reset 6527029 --hard \
    && cd ../ \
    && mkdir GenOO \
    && cp -r GenOO_git/lib/GenOO/* GenOO/

# Install specific version of Perl module https://stackoverflow.com/questions/260593/how-can-i-install-a-specific-version-of-a-set-of-perl-modules
RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm --force CLIPSeqTools@0.1.9  \
    && cp -r /usr/local/TERA-Seq_manuscript/misc/GenOOx/* perl-virtualenv/teraseq/lib/perl5/GenOOx/

####################################################################################################
# ### Nanopolish
# # Default version
# RUN git clone --recursive https://github.com/jts/nanopolish.git \
#     && mv nanopolish nanopolish-480fc85 \
#     && cd nanopolish-480fc85/ \
#     && git reset 480fc85 --hard \
#     && sed -i 's#http://bitbucket.org/eigen/eigen/get/$(EIGEN_VERSION).tar.bz2#https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.bz2#' Makefile \
#     && sed -i 's/tar -xjf $(EIGEN_VERSION).tar.bz2/tar -xjf eigen-$(EIGEN_VERSION).tar.bz2/' Makefile \
#     && sed -i 's/eigen-eigen-\*/eigen-$(EIGEN_VERSION)/' Makefile \
#     && rm -rf fast5 \
#     && git clone https://github.com/mateidavid/fast5.git \
#     && cd fast5/ \
#     && git reset 18d6e34 --hard \
#     && cd ../ \
#     && rm -rf htslib \
#     && git clone --recursive https://github.com/samtools/htslib.git \
#     && cd htslib/ \
#     && git reset 3dc96c5 --hard \
#     && cd ../ \
#     && make \
#     && ln -s $(pwd)/nanopolish /usr/local/miniconda3/envs/teraseq/bin/nanopolish

# # New version with polya hmm scripts
# RUN git clone --recursive https://github.com/jts/nanopolish.git \
#     && mv nanopolish nanopolish-ab9722b \
#     && cd nanopolish-ab9722b/ \
#     && git reset ab9722b --hard
####################################################################################################

### Other dependencies
# Make sure to activate Conda
SHELL ["conda", "run", "-n", "teraseq", "/bin/bash", "-c"]

## GeneCycle
# TODO: Install specific version of R package https://support.posit.co/hc/en-us/articles/219949047-Installing-older-versions-of-packages
RUN Rscript -e 'install.packages("GeneCycle", repos="https://cloud.r-project.org")'

## Cutadapt
RUN mkdir cutadapt-2.5 \
    && cd cutadapt-2.5/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install cutadapt==2.5 pysam numpy pandas matplotlib seaborn \
    && which cutadapt

## DeepTools
RUN mkdir deepTools-3.5.0 \
    && cd deepTools-3.5.0/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install wheel \
    && pip3 install deeptools==3.5.0 \
    && deeptools --version

## ONT-Fast5-API
RUN mkdir ont-fast5-api \
    && cd ont-fast5-api/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && pip install ont-fast5-api==3.3.0 h5py seaborn

## Jvarkit
RUN git clone "https://github.com/lindenb/jvarkit.git" \
    && mv jvarkit jvarkit-014d3e9 \
    && cd jvarkit-014d3e9/ \
    && git reset 014d3e9 --hard \
    && ./gradlew biostar84452 \
    && mkdir $CONDA_PREFIX/share/jvarkit \
    && ln -s $(pwd)/dist/biostar84452.jar /usr/local/miniconda3/envs/teraseq/share/jvarkit/remove-softlip.jar

# Set permissions the same for all users to avoid Singularity "issue" with running as root (Docker) vs running as user who executed the command (Singularity)
RUN chmod -R a=u /usr/local/TERA-Seq_manuscript/tools

### Add utils dir to PATH
ENV PATH "/usr/local/TERA-Seq_manuscript/tools/utils:${PATH}"

WORKDIR /usr/local/TERA-Seq_manuscript
