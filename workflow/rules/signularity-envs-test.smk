rule conda_test:
    output:
        "results/singularity-conda-test.txt"
    shell:
        '''
        {activ_conda}
        echo $CONDA_PREFIX > {output}
        samtools --version >> {output}
        '''

rule perl_test:
    output:
        "results/singularity-perl-test.txt"
    shell:
        '''
        {activ_perl}
        which perl > {output}
        /usr/local/TERA-Seq_manuscript/tools/utils/fastq-sanitize-header -h >> {output}
        '''

rule cutadapt_test:
    output:
        "results/singularity-cutadapt-test.txt"
    shell:
        '''
        {activ_cutadapt}
        which python > {output}
        cutadapt --version >> {output}
        '''

rule deeptools_test:
    output:
        "results/singularity-deeptools-test.txt"
    shell:
        '''
        {activ_deeptools}
        which python > {output}
        computeMatrix --version >> {output}
        '''