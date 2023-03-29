rule map_trans_polya_minimap2:
    input:
        fastq="data/samples/{sample}/fastq/reads.1.sanitize.noribo.fastq.gz",
        mmi_trans=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['mmi_trans'],
    output:
        bam=["data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam",
        "data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam"],
    params:
        k=12,
        secondary="yes",
        sec_to_prim=1.0,
        mem_mb=768, # samtools sort sets -m per thread
    threads: 32
    shell:
        '''
        {activ_conda}

        minimap2 \
            -a \
            -x map-ont \
            -k {params.k} \
            -u f \
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

        ln -sf $(basename {output.bam[0]}) \
            {output.bam[1]}
        '''

rule map_genome_minimap2:
    input:
#        fastq="data/samples/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        fastq="data/samples/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
        mmi_genome=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['mmi_genome'],
    output:
        bam="data/samples/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    params:
        k=12,
        secondary="yes",
        sec_to_prim=1.0,
        mem_mb=768, # samtools sort sets -m per thread
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
        bam="data/samples/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",    
    output:
        bai="data/samples/{sample}/align/reads.1.sanitize.toGenome.sorted.bam.bai",    
    threads: 1
    shell:
        '''
        {activ_conda}

        samtools index -@ {threads} \
            {input.bam}
        '''

rule index_bam_trans_polya:
    input:
        bam="data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam",
    output:
        bai=["data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam.bai",
        "data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam.bai"],
    threads: 1
    shell:
        '''
        {activ_conda}

        samtools index -@ {threads} \
            {input.bam} \
        
        ln -sf $(basename {output[0]}) \
            {output[1]}
        '''