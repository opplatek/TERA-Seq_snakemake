rule conda_test:
    output:
        "results/singularity-conda-test.txt"
    shell:
        '''
        source /root/miniconda3/bin/activate
        conda activate teraseq
        echo $CONDA_PREFIX > {output}
        samtools --version >> {output}
        '''

rule perl_test:
    output:
        "results/singularity-perl-test.txt"
    shell:
        '''
        set +eu
        source /root/TERA-Seq_manuscript/tools/perl-virtualenv/teraseq/bin/activate
        set -eu
        which perl > {output}
        /root/TERA-Seq_manuscript/tools/utils/fastq-sanitize-header -h >> {output}
        '''

rule cutadapt_test:
    output:
        "results/singularity-cutadapt-test.txt"
    shell:
        '''
        source /root/miniconda3/bin/activate
        conda activate teraseq
        source /root/TERA-Seq_manuscript/tools/cutadapt-2.5/venv/bin/activate
        which python > {output}
        cutadapt --version >> {output}
        '''

rule deeptools_test:
    output:
        "results/singularity-deeptools-test.txt"
    shell:
        '''
        source /root/miniconda3/bin/activate
        conda activate teraseq
        source /root/TERA-Seq_manuscript/tools/deepTools-3.5.0/venv/bin/activate
        which python > {output}
        computeMatrix --version >> {output}
        '''