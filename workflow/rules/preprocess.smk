def get_adapter(libtypes, adapters, sample):
    adapter=adapters[libtypes[sample]]
    return adapter

rule adapter_remove_rel5:
    input:
        "data/samples/{sample}/fastq/reads.1.sanitize.fastq.gz"
    output:
        wadapt="data/samples/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
        woadapt="data/samples/{sample}/fastq/reads.1.sanitize.wo_rel5.fastq.gz",
    params:
#        rel5long="XAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        rel5long=lambda wildcards: get_adapter(LIBTYPES, ADAPTERS, wildcards.sample),        
        overlap=31,
        minlen=25,
        errorrate=0.29,
    log:
        "data/samples/{sample}/log/cutadapt.log"
    shell:
        '''
        {activ_cutadapt}

        cutadapt \
            -g {params.rel5long} \
            --overlap {params.overlap} \
            --minimum-length {params.minlen} \
            --error-rate {params.errorrate} \
            --output {output.wadapt} \
            --untrimmed-output {output.woadapt} \
            {input} \
            &> {log}
        '''

rule fastq_adapter_merge_rel5:
    input: 
        wadapt="data/samples/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
        woadapt="data/samples/{sample}/fastq/reads.1.sanitize.wo_rel5.fastq.gz",
    output:
        wadapt_names="data/samples/{sample}/fastq/reads.1.sanitize.w_rel5.names.txt",
        woadapt_names="data/samples/{sample}/fastq/reads.1.sanitize.wo_rel5.names.txt",
        merged="data/samples/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",    
    shell:
        '''
        zcat {input.wadapt} \
            | paste - - - - | cut -f1 | sed 's/^@//g' \
            > {output.wadapt_names}
        zcat {input.woadapt} \
            | paste - - - - | cut -f1 | sed 's/^@//g' \
            > {output.woadapt_names}

        cat {input.wadapt} {input.woadapt} \
            > {output.merged}
        '''

rule ribosomal_map_minimap2:
    input:
        fastq="data/samples/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        mmi_wribo=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['mmi_wribo']
#        mmi_wribo="data/hg38/minimap2.17/ensembl-transcripts-wRibo.k12.mmi",
    output:
        temp("data/samples/{sample}/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam"),
    params:
        k=12,
        secondary="yes",
        mem_mb=768, # samtools sort sets -m per thread
    threads: 32
    shell:
        '''
        {activ_conda}

        echo {input.mmi_wribo}

        minimap2 \
            -a \
            -x map-ont \
            -k {params.k} \
            -p 1 \
            -u f \
            -t {threads} \
            --secondary={params.secondary} \
            {input.mmi_wribo} \
            {input.fastq} \
        | samtools view -b - \
        | samtools sort -@ {threads} -m {params.mem_mb}M - \
            > {output}   
        '''

rule ribosomal_extract:
    input:
        "data/samples/{sample}/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam",
    output:
        sam=temp("data/samples/{sample}/align/reads.1.sanitize.toRibosomal.sorted.sam"),
        ribo_names=temp("data/samples/{sample}/align/reads.1.sanitize.toRibosomal.sorted.reads.txt"),  
        bam="data/samples/{sample}/align/reads.1.sanitize.toRibosomal.sorted.bam",      
    threads: 32
    shell:
        '''
        {activ_conda}

        samtools view -H {input} \
            > {output.sam}

        samtools view -@ {threads} -F4 {input} \
            | grep -v -P "\\tENST" \
            >> {output.sam}

        cat {output.sam} \
            | cut -f1 | awk '!x[$0]++' \
            | grep -P -v "^@HD|^@SQ|^@PG" \
            > {output.ribo_names}

        samtools view -@ {threads} -bh {output.sam} | samtools sort -@ {threads} - \
            > {output.bam}
        '''

rule ribosomal_remove:
    input:
        ribo_names="data/samples/{sample}/align/reads.1.sanitize.toRibosomal.sorted.reads.txt",
        fastq="data/samples/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
    output:
        fastq="data/samples/{sample}/fastq/reads.1.sanitize.noribo.fastq.gz",
    shell:
        '''
        {activ_conda}

        seqkit grep -nvf {input.ribo_names} \
            {input.fastq} \
            -o {output.fastq}
        '''