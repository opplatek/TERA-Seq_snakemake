# IMPORTANT: Sql isnt't compatible with Snakemake because we are constantly updating a file. Snakemake then thinks the input has changed and will try to redo all the rules. We would have to omit sql db file from the inputs.

### Transcriptome
rule sam_to_sqlite_trans:
    input:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        db=samplesdir + "/{sample}/db/sqlite.db",
        done=samplesdir + "/{sample}/log/sam_to_sqlite_trans.done",
    params:
        filter_flags="-F 4 -F 16 -F 2048",
        table="transcr",
        conda_path=CONDA_PATH,
    shell:
        '''
        {activ_perl}

        {params.conda_path}/samtools view -h {params.filter_flags} {input.bam} \
            | sam_to_sqlite \
                --database {output.db} \
                --table {params.table} \
                --records_class GenOOx::Data::File::SAMminimap2::Record \
                --drop \
            && touch {output.done}
        '''


rule annotate_sqlite_trans_mrna:
    input:
        done=samplesdir + "/{sample}/log/sam_to_sqlite_trans.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['bed_mrna'],
    output:
        done=samplesdir + "/{sample}/log/annot_trans_mrna.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="coding_transcript",
    shell:
        '''
        {activ_perl}

        clipseqtools-preprocess annotate_with_file \
            --database {params.db} \
            --table {params.table} \
            --a_file {input.bed} \
            --column {params.column} \
        && touch {output.done}
        '''


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_ncrna with:
    input:
        done=samplesdir + "/{sample}/log/annot_trans_mrna.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['bed_ncrna'],
    output:
        done=samplesdir + "/{sample}/log/annot_trans_ncrna.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="noncoding_transcript",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr5 with:
    input:
        done=samplesdir + "/{sample}/log/annot_trans_ncrna.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['bed_utr5'],
    output:
        done=samplesdir + "/{sample}/log/annot_trans_utr5.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="utr5",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_cds with:
    input:
        done=samplesdir + "/{sample}/log/annot_trans_utr5.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['bed_cds'],
    output:
        done=samplesdir + "/{sample}/log/annot_trans_cds.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="cds",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr3 with:
    input:
        done=samplesdir + "/{sample}/log/annot_trans_cds.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['bed_utr3'],
    output:
        done=samplesdir + "/{sample}/log/annot_trans_utr3.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="utr3",


#rule annotate_sqlite_trans_adapter_rel5:
rule annotate_sqlite_trans_adapter:
    input:
        done=samplesdir + "/{sample}/log/annot_trans_utr3.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
#        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
    output:
#        done=samplesdir + "/{sample}/log/annot_trans_rel5.done",
        done=samplesdir + "/{sample}/log/annot_trans_adapt.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
#        column="rel5",
        column=lambda wildcards: get_adaptside(LIBTYPES, wildcards.sample),
        column_bind="qname",
    shell:
        '''
        {activ_conda}

        annotate-sqlite-with-fastq \
            --database {params.db} \
            --db_col_bind {params.column_bind} \
            --db_col_add {params.column} \
            --db_tables {params.table} \
            --ifile {input.wadapt} \
        && touch {output.done}
        '''


### Genome
rule sam_to_sqlite_genome:
    input:
#        done=samplesdir + "/{sample}/log/annot_trans_rel5.done",
        done=samplesdir + "/{sample}/log/annot_trans_adapt.done",
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        done=samplesdir + "/{sample}/log/sam_to_sqlite_genome.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        filter_flags="-F 4 -F 2048",
        table="genome",
        conda_path=CONDA_PATH,
    shell:
        '''
        {activ_perl}

        {params.conda_path}/samtools view -h {params.filter_flags} {input.bam} \
            | sam_to_sqlite \
                --database {params.db} \
                --table {params.table} \
                --records_class GenOOx::Data::File::SAMminimap2::Record \
                --drop \
            && touch {output.done}
        '''


rule annotate_sqlite_genome_gtf_polya:
    input:
        done=samplesdir + "/{sample}/log/sam_to_sqlite_genome.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        gtf=lambda wildcards: get_refs(datadir, ASSEMBLIES, wildcards.sample)['gtf_polya'],
    output:
        done=samplesdir + "/{sample}/log/annot_genome_gtf.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="genome",
    shell:
        '''
        {activ_perl}

        clipseqtools-preprocess annotate_with_genic_elements \
            --database {params.db} \
            --table {params.table} \
            --gtf {input.gtf} \
        && touch {output.done}
        '''


#use rule annotate_sqlite_trans_adapter_rel5 as annotate_sqlite_genome_adapter_rel5 with:
use rule annotate_sqlite_trans_adapter as annotate_sqlite_genome_adapter with:
    input:
        done=samplesdir + "/{sample}/log/annot_genome_gtf.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
#        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",
        wadapt=samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.fastq.gz",
    output:
#        done=samplesdir + "/{sample}/log/annot_genome_rel5.done",
        done=samplesdir + "/{sample}/log/annot_genome_adapt.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="genome",
#        column="rel5",
        column=lambda wildcards: get_adaptside(LIBTYPES, wildcards.sample),
        column_bind="qname",

