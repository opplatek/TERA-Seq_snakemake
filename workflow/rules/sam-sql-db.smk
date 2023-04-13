# IMPORTANT: Sql isnt't compatible with Snakemake because we are constantly updating a file. Snakemake then thinks the input has changed and will try to redo all the rules. We would have to omit sql db file from the inputs.

### Transcriptome
rule sam_to_sqlite_trans:
    input:
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        db=samplesdir + "/{sample}/db/sqlite.db",
        done=samplesdir + "/{sample}/db/sam_to_sqlite_trans.done",
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
        done=samplesdir + "/{sample}/db/sam_to_sqlite_trans.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['bed_mrna'],
    output:
        done=samplesdir + "/{sample}/db/annot_trans_mrna.done",
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
        done=samplesdir + "/{sample}/db/annot_trans_mrna.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['bed_ncrna'],
    output:
        done=samplesdir + "/{sample}/db/annot_trans_ncrna.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="noncoding_transcript",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr5 with:
    input:
        done=samplesdir + "/{sample}/db/annot_trans_ncrna.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['bed_utr5'],
    output:
        done=samplesdir + "/{sample}/db/annot_trans_utr5.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="utr5",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_cds with:
    input:
        done=samplesdir + "/{sample}/db/annot_trans_utr5.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['bed_cds'],
    output:
        done=samplesdir + "/{sample}/db/annot_trans_cds.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="cds",


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr3 with:
    input:
        done=samplesdir + "/{sample}/db/annot_trans_cds.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['bed_utr3'],
    output:
        done=samplesdir + "/{sample}/db/annot_trans_utr3.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
        column="utr3",


rule annotate_sqlite_trans_5tera:
    input:
        done=samplesdir + "/{sample}/db/annot_trans_utr3.done",
        wadapt = samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.5tera.fastq.gz",      
    output:
        done=samplesdir + "/{sample}/db/annot_trans_adapt.5tera.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
#        column=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample), # RSQLite doesn't like column names starting with a nubmer
        column="rel5",
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


use rule annotate_sqlite_trans_5tera as annotate_sqlite_trans_tera3 with:
    input:
        done=lambda wildcards: get_annot_trans_tera3_input(LIBTYPES, wildcards.sample),
        wadapt = samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.tera3.fastq.gz",      
    output:
        done=samplesdir + "/{sample}/db/annot_trans_adapt.tera3.done",
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="transcr",
#        column=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample), # RSQLite doesn't like column names starting with a nubmer
        column="rel3",
        column_bind="qname",

### Genome
rule sam_to_sqlite_genome:
    input:
#        done=samplesdir + "/{sample}/db/annot_trans_adapt.done",
        done=lambda wildcards: get_sam_to_sqlite_genome_input(LIBTYPES, wildcards.sample),
        bam=samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        done=samplesdir + "/{sample}/db/sam_to_sqlite_genome.done",
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


rule annotate_sqlite_genome_gtf:
    input:
        done=samplesdir + "/{sample}/db/sam_to_sqlite_genome.done",
        # db=samplesdir + "/{sample}/db/sqlite.db",
        gtf=lambda wildcards: get_refs_trans(datadir, ASSEMBLIES, PROTOCOLS, wildcards.sample)['gtf'],
    output:
        done=samplesdir + "/{sample}/db/annot_genome_gtf.done",
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
use rule annotate_sqlite_trans_5tera as annotate_sqlite_genome_5tera with:
    input:
        done=samplesdir + "/{sample}/db/annot_genome_gtf.done",
#        wadapt = lambda wildcards: get_trim_merge(LIBTYPES, wildcards.sample)['wadapt']
        wadapt = samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.5tera.fastq.gz",      
    output:
        done=samplesdir + "/{sample}/db/annot_genome_adapt.5tera.done",        
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="genome",
#        column=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample), # RSQLite doesn't like column names starting with a nubmer
        column="rel5",
        column_bind="qname",


use rule annotate_sqlite_genome_5tera as annotate_sqlite_genome_tera3 with:
    input:
#        done=samplesdir + "/{sample}/db/annot_genome_gtf.done",
        done=lambda wildcards: get_annot_genome_tera3_input(LIBTYPES, wildcards.sample),
#        wadapt = lambda wildcards: get_trim_merge(LIBTYPES, wildcards.sample)['wadapt']
        wadapt = samplesdir + "/{sample}/fastq/reads.1.sanitize.w_adapt.tera3.fastq.gz",      
    output:
        done=samplesdir + "/{sample}/db/annot_genome_adapt.tera3.done",        
    params:
        db=samplesdir + "/{sample}/db/sqlite.db",
        table="genome",
#        column=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample), # RSQLite doesn't like column names starting with a nubmer
        column="rel3",
        column_bind="qname",


rule sql_annotate:
    input:
        lambda wildcards: get_sqldb_done(LIBTYPES, wildcards.sample)
    output:
        samplesdir + "/{sample}/db/annot_sqldb.done",
