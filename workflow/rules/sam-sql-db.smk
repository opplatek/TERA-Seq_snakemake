ruleorder: sam_to_sqlite_trans > sam_to_sqlite_genome

### Transcriptome
rule sam_to_sqlite_trans:
    input:
        bam="data/samples/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        db="data/samples/{sample}/db/sqlite.db",
        done="data/samples/{sample}/db/sam_to_sqlite_trans.done",
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
        done="data/samples/{sample}/db/sam_to_sqlite_trans.done",
        db="data/samples/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['bed_mrna'], 
    output:
        done="data/samples/{sample}/db/annot_trans_mrna.done",
    params:
        table="transcr",  
        column="coding_transcript",   
    shell:
        '''
        {activ_perl}

        clipseqtools-preprocess annotate_with_file \
            --database {input.db} \
            --table {params.table} \
            --a_file {input.bed} \
            --column {params.column} \
        && touch {output.done}
        '''

use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_ncrna with:
    input:
        done="data/samples/{sample}/db/annot_trans_mrna.done",
        db="data/samples/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['bed_ncrna'], 
    output:
        done="data/samples/{sample}/db/annot_trans_ncrna.done",
    params:
        table="transcr",
        column="noncoding_transcript",   


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr5 with:
    input:
        done="data/samples/{sample}/db/annot_trans_ncrna.done",
        db="data/samples/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['bed_utr5'],
    output:
        done="data/samples/{sample}/db/annot_trans_utr5.done",
    params:
        table="transcr",
        column="utr5", 


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_cds with:
    input:
        done="data/samples/{sample}/db/annot_trans_utr5.done",
        db="data/samples/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['bed_cds'],  
    output:
        done="data/samples/{sample}/db/annot_trans_cds.done",
    params:
        table="transcr",    
        column="cds",     


use rule annotate_sqlite_trans_mrna as annotate_sqlite_trans_utr3 with:
    input:
        done="data/samples/{sample}/db/annot_trans_cds.done",
        db="data/samples/{sample}/db/sqlite.db",
        bed=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['bed_utr3'],    
    output:
        done="data/samples/{sample}/db/annot_trans_utr3.done",
    params:
        table="transcr",
        column="utr3",         


rule annotate_sqlite_trans_adapter_rel5:
    input:
        done="data/samples/{sample}/db/annot_trans_utr3.done",
        db="data/samples/{sample}/db/sqlite.db",
        wadapt="data/samples/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",            
    output:
        done="data/samples/{sample}/db/annot_trans_rel5.done",    
    params:
        table="transcr",
        column="rel5",
        column_bind="qname",
    shell:
        '''
        {activ_conda}

        annotate-sqlite-with-fastq \
            --database {input.db} \
            --db_col_bind {params.column_bind} \
            --db_col_add {params.column} \
            --db_tables {params.table} \
            --ifile {input.wadapt} \
        && touch {output.done}
        '''


### Genome
use rule sam_to_sqlite_trans as sam_to_sqlite_genome with:
    input:
        done="data/samples/{sample}/db/annot_trans_rel5.done",
        bam="data/samples/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        db="data/samples/{sample}/db/sqlite.db",
        done="data/samples/{sample}/db/sam_to_sqlite_genome.done",
    params:
        filter_flags="-F 4 -F 2048",
        table="genome",
        conda_path=CONDA_PATH,


rule annotate_sqlite_genome_gtf_polya:
    input:
        done="data/samples/{sample}/db/sam_to_sqlite_genome.done",
        db="data/samples/{sample}/db/sqlite.db",
        gtf=lambda wildcards: get_refs(ASSEMBLIES, wildcards.sample)['gtf_polya'],
    output:
        done="data/samples/{sample}/db/annot_genome_gtf.done",
    params:
        table="genome",
    shell:  
        '''      
        {activ_perl}

        clipseqtools-preprocess annotate_with_genic_elements \
            --database {input.db} \
            --table {params.table} \
            --gtf {input.gtf} \
        && touch {output.done}
        '''

use rule annotate_sqlite_trans_adapter_rel5 as annotate_sqlite_genome_adapter_rel5 with:
    input:
        done="data/samples/{sample}/db/annot_genome_gtf.done",
        db="data/samples/{sample}/db/sqlite.db",
        wadapt="data/samples/{sample}/fastq/reads.1.sanitize.w_rel5.fastq.gz",            
    output:
        done="data/samples/{sample}/db/annot_genome_rel5.done",
    params:
        table="genome",
        column="rel5",
        column_bind="qname",        

