rule fastq_sanitize:
    input:
        "data/samples/{sample}/fastq/reads.1.fastq.gz"
    output:
        "data/samples/{sample}/fastq/reads.1.sanitize.fastq.gz"
    shell:
        '''
        {activ_perl}

        zcat {input} \
            | fastq-sanitize-header --input - --delim : --keep 0 \
            | gzip -c \
            > {output}
        '''
