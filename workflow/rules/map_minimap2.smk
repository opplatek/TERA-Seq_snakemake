rule map_trans_polya_minimap2:
    input:
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.noribo.fastq.gz",
        mmi_trans=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['mmi_trans'],
    output:
        bam=[samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam",
             samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam"],
    params:
        k=12,
        secondary="yes",
        sec_to_prim=1.0,
        mem_mb=samtools_sort_mb # samtools sort sets -m per thread; 768 is default
    resources:
        map_jobs=4,
        mem_mb=get_mem_mb
    threads: 32
    shell:
        '''
        {activ_conda}

        minimap2 \
            -a \
            -x map-ont \
            -k {params.k} \
            --for-only \
            -t {threads} \
            --secondary={params.secondary} \
            -p {params.sec_to_prim} \
            {input.mmi_trans} \
            {input.fastq} \
        | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
        | grep -v "SA:Z:" \
        | sam-count-secondary --tag X0 \
        | samtools view -b - \
        | samtools sort -@ {threads} -m {params.mem_mb}M - \
        > {output.bam[0]}

        # In the TERA-Seq paper, we used '-u f'. However, this doesn't mean the read alignment would be 'forced' to align to the forward strand (TERA-Seq is forward-strand specific)
        # '--for-only' is the correct parameter for Minimap2 to consider only forward-strand alignments

        ln -sf $(basename {output.bam[0]}) \
            {output.bam[1]}
        '''

rule map_genome_minimap2:
    input:
#        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
        mmi_genome=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['mmi_genome'],
    output:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    params:
        k=12,
        secondary="yes",
        sec_to_prim=1.0,
        mem_mb=samtools_sort_mb # samtools sort sets -m per thread; 768 is default
    resources:
        map_jobs=4,
        mem_mb=get_mem_mb
    threads: 32
    shell:
        '''
        {activ_conda}

        minimap2 \
            -a \
            -x splice \
            -k {params.k} \
            -u b \
            -t {threads} \
            --secondary={params.secondary} \
            -p {params.sec_to_prim} \
            {input.mmi_genome} \
            {input.fastq} \
        | add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
        | grep -v "SA:Z:" \
        | sam-count-secondary --tag X0 \
        | samtools view -b - \
        | samtools sort -@ {threads} -m {params.mem_mb}M - \
        > {output.bam}
        '''

rule index_bam_genome:
    input:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        bai=samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam.bai",
    threads: 1
    shell:
        '''
        {activ_conda}

        samtools index -@ {threads} \
            {input.bam}
        '''

rule index_bam_trans_polya:
    input:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam",
    output:
        bai=[samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam.bai",
             samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam.bai"],
    threads: 1
    shell:
        '''
        {activ_conda}

        samtools index -@ {threads} \
            {input.bam} \

        ln -sf $(basename {output[0]}) \
            {output[1]}
        '''
