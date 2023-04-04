#rule adapter_remove_rel5:
rule adapter_remove:
    input:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.fastq.gz"
    output:
#        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
#        woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_rel5.fastq.gz",
        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
        woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.fastq.gz",
    params:
#        adapter="XAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        adapter=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['sequence'],
        libtype=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample),
        # side="-g"
        # overlap=31,
        # minlen=25,
        # errorrate=0.29,
        side=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['side'],
        overlap=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['overlap'],
        minlen=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['minlen'],
        errorrate=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['errorrate'],
    log:
        samplesdir + "/{sample}/log/cutadapt.log"
    shell:
        '''
        {activ_cutadapt}

        cutadapt \
            {params.side} {params.adapter} \
            --overlap {params.overlap} \
            --minimum-length {params.minlen} \
            --error-rate {params.errorrate} \
            --output {output.wadapt} \
            --untrimmed-output {output.woadapt} \
            {input} \
            &> {log}
        '''


#rule fastq_adapter_names_rel5:
rule fastq_adapter_names:
    input:
#       wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
#       woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_rel5.fastq.gz",
       wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
       woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.fastq.gz",
    output:
        # wadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.names.txt",
        # woadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_rel5.names.txt",
        wadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.names.txt",
        woadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.names.txt",
    shell:
        '''
        zcat {input.wadapt} \
            | paste - - - - | cut -f1 | sed 's/^@//g' \
            > {output.wadapt_names}
        zcat {input.woadapt} \
            | paste - - - - | cut -f1 | sed 's/^@//g' \
            > {output.woadapt_names}
        '''

#rule fastq_adapter_merge_rel5:
rule fastq_adapter_merge:
    input:
    #    wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
    #    woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_rel5.fastq.gz",
       wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
       woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.fastq.gz",
    output:
#        merged=samplesdir + "/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        merged=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
    shell:
        '''
        cat {input.wadapt} {input.woadapt} \
            > {output.merged}
        '''


rule ribosomal_map_minimap2:
    input:
#        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
#        mmi_wribo=datadir + "/hg38/minimap2.17/ensembl-transcripts-wRibo.k12.mmi",
        mmi_wribo=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['mmi_wribo']
    output:
        temp(samplesdir + "/{sample}/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam"),
    params:
        k=12,
        secondary="yes",
        mem_mb=samtools_sort_mb # samtools sort sets -m per thread; 768 is default
    resources:
        mem_mb=get_mem_mb
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
        samplesdir + "/{sample}/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam",
    output:
        sam=temp(samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.sam"),
        ribo_names=temp(samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.reads.txt"),
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.bam",
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
        ribo_names=samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.reads.txt",
#        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.rel5_trim.fastq.gz",
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
    output:
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.noribo.fastq.gz",
    shell:
        '''
        {activ_conda}

        seqkit grep -nvf {input.ribo_names} \
            {input.fastq} \
            -o {output.fastq}
        '''
