rule adapter_remove_5tera:
    input:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.fastq.gz"
    output:
        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.5tera.fastq.gz",
        woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.5tera.fastq.gz",
        merged=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.5tera.fastq.gz",
        wadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.5tera.names.txt",
        woadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.5tera.names.txt",
        trimlog=samplesdir + "/{sample}/log/cutadapt.5tera.log",
    params:
        # adapter=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['sequence'],
        # side=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['side'],
        # overlap=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['overlap'],
        # minlen=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['minlen'],
        # errorrate=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['errorrate'],
        adapter=lambda wildcards: get_cutadapt_settings(ADAPTERS, '5tera')['sequence'],
        side=lambda wildcards: get_cutadapt_settings(ADAPTERS, '5tera')['side'],
        overlap=lambda wildcards: get_cutadapt_settings(ADAPTERS, '5tera')['overlap'],
        minlen=lambda wildcards: get_cutadapt_settings(ADAPTERS, '5tera')['minlen'],
        errorrate=lambda wildcards: get_cutadapt_settings(ADAPTERS, '5tera')['errorrate'],
#    log:
#        samplesdir + "/{sample}/log/cutadapt.5tera.log",
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
                &> {output.trimlog}

            zcat {output.wadapt} \
                | paste - - - - | cut -f1 | sed 's/^@//g' \
                > {output.wadapt_names}
            zcat {output.woadapt} \
                | paste - - - - | cut -f1 | sed 's/^@//g' \
                > {output.woadapt_names}

            cat {output.wadapt} {output.woadapt} \
                > {output.merged}                
        '''


rule adapter_remove_tera3:
    input:
#        samplesdir + "/{sample}/fastq/reads.1.sanitize.fastq.gz"
        lambda wildcards: get_tera3_remove_input(LIBTYPES, wildcards.sample)
    output:
        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.tera3.fastq.gz",
        woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.tera3.fastq.gz",
        merged=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.tera3.fastq.gz",
        wadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.tera3.names.txt",
        woadapt_names=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.tera3.names.txt",
        fastq_temp_last = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_last.tera3.fastq.gz"),
        fastq_temp_short = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_short.tera3.fastq.gz"),
        fastq_temp_merged = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_merged.tera3.fastq.gz"),
        fastq_temp_merged_w_tera3 = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_merged.w_tera3.fastq.gz"),
        fastq_temp_merged_wo_tera3 = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_merged.wo_tera3.fastq.gz"),
        fastq_temp_merged_empty = temp(samplesdir + "/{sample}/fastq/reads.1.sanitize.tmp_merged.ishouldbeempty.tera3.fastq.gz"),
        trimlog=samplesdir + "/{sample}/log/cutadapt.tera3.log",
        trimlog_fulllen=samplesdir + "/{sample}/log/cutadapt.tera3.full-length.log",        
    params:
        # adapter=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['sequence'],
        # side=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['side'],
        # overlap=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['overlap'],
        # minlen=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['minlen'],
        # errorrate=lambda wildcards: get_cutadapt_settings(LIBTYPES, ADAPTERS, wildcards.sample)['errorrate'],
        adapter=lambda wildcards: get_cutadapt_settings(ADAPTERS, 'tera3')['sequence'],
        side=lambda wildcards: get_cutadapt_settings(ADAPTERS, 'tera3')['side'],
        overlap=lambda wildcards: get_cutadapt_settings(ADAPTERS, 'tera3')['overlap'],
        minlen=lambda wildcards: get_cutadapt_settings(ADAPTERS, 'tera3')['minlen'],
        errorrate=lambda wildcards: get_cutadapt_settings(ADAPTERS, 'tera3')['errorrate'],        
        tera3_nucl_search=200, # How many nucleotides from the 3' end of TERA3 libraries to scan for adapter
    # log:
    #     samplesdir + "/{sample}/log/cutadapt.tera3.log",
    #     samplesdir + "/{sample}/log/cutadapt.tera3.full-length.log",
    shell:
        '''
        {activ_cutadapt}

            # Get only last N nucleotides to scan for adapter
            gunzip -c {input} \
                | seqkit seq -m $[{params.tera3_nucl_search}+1] \
                | seqkit subseq -r -{params.tera3_nucl_search}:-1 \
                | gzip -c > {output.fastq_temp_last}

            # Get shorter reads just not to exclude them at the start
            gunzip -c {input} \
                | seqkit seq -M {params.tera3_nucl_search} \
                | gzip -c > {output.fastq_temp_short}

            cat {output.fastq_temp_last} {output.fastq_temp_short} \
                > {output.fastq_temp_merged}

            # First round to identify reads with adapters in the last N nucleotides
            cutadapt \
                -a GTGTCAGTCACTTCCA \
                --overlap {params.overlap} \
                --minimum-length 0 \
                --error-rate {params.errorrate} \
                --untrimmed-output {output.fastq_temp_merged_wo_tera3} \
                --output {output.fastq_temp_merged_w_tera3} \
                {output.fastq_temp_merged} \
                &> {output.trimlog}

            zcat {output.fastq_temp_merged_w_tera3} \
                | paste - - - - | cut -f1 | sed 's/^@//g' \
                > {output.wadapt_names}
            zcat {output.fastq_temp_merged_wo_tera3} \
                | paste - - - - | cut -f1 | sed 's/^@//g' \
                > {output.woadapt_names}

            zcat {input} \
                | seqtk subseq - {output.wadapt_names} \
                | gzip -c > {output.fastq_temp_merged_w_tera3}
            zcat {input} \
                | seqtk subseq - {output.woadapt_names} \
                | seqkit seq -m {params.minlen} | gzip -c > {output.woadapt}

            cutadapt \
                -a GTGTCAGTCACTTCCA \
                --overlap {params.overlap} \
                --minimum-length {params.minlen} \
                --error-rate {params.errorrate} \
                --untrimmed-output {output.fastq_temp_merged_empty} \
                --output {output.wadapt} \
                {output.fastq_temp_merged_w_tera3} \
                &> {output.trimlog_fulllen}
            
            cat {output.wadapt} {output.woadapt} \
                > {output.merged}
        '''


rule fastq_adapter_remove_sync:
    input:
        lambda wildcards: get_adapter_remove_sync_input(LIBTYPES, wildcards.sample)
    output:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
    shell:
        '''
        ln -sf $(basename {input}) {output}
        '''


# rule fastq_adapter_merge:
#     input:
# #        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
# #        woadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.wo_adapt.fastq.gz",
#         unpack(lambda wildcards: get_trim_merge(LIBTYPES, wildcards.sample)),
#     output:
#         merged=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
#     shell:
#         '''
#         cat {input.wadapt} {input.woadapt} \
#             > {output.merged}
#         '''


rule ribosomal_map_minimap2:
    input:
        fastq=samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
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

        minimap2 \
            -a \
            -x map-ont \
            -k {params.k} \
            -p 1 \
            --for-only \
            -t {threads} \
            --secondary={params.secondary} \
            {input.mmi_wribo} \
            {input.fastq} \
        | samtools view -b - \
        | samtools sort -@ {threads} -m {params.mem_mb}M - \
            > {output}

        # In the TERA-Seq paper, we used '-u f'. However, this doesn't mean the read alignment would be 'forced' to align to the forward strand (TERA-Seq is forward-strand specific)
        # '--for-only' is the correct parameter for Minimap2 to consider only forward-strand alignmen
        '''


rule ribosomal_extract:
    input:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toEnsembl-transcripts-wRibo.sorted.bam",
        names_rrna=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['names_rrna'],
    output:
        sam=temp(samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.sam"),
        ribo_names=temp(samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.reads.txt"),
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toRibosomal.sorted.bam",
    threads: 32
    shell:
        '''
        {activ_conda}

        samtools view -H {input.bam} \
            > {output.sam}

        samtools view -@ {threads} -F 4 {input.bam} \
            | grep -wFf {input.names_rrna} \
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
