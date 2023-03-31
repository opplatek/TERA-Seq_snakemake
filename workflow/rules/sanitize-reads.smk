rule fastq_sanitize:
    input:
        samplesdir + "/{sample}/fastq/reads.1.fastq.gz"
    output:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.fastq.gz"
    shell:
        '''
        {activ_perl}

        zcat {input} \
            | fastq-sanitize-header --input - --delim : --keep 0 \
            | gzip -c \
            > {output}
        '''
