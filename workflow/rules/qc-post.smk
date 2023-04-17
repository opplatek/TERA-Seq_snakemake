# rule cutadapt_parse_lens:
#     input:
#         samplesdir + "/{sample}/log/cutadapt.log",
#     output:
#         samplesdir + "/{sample}/log/cutadapt.len.tsv",
#     shell:
#         '''
#         {activ_conda}
#
#         python3 workflow/scripts/parse-cutadapt-lens.py {input} {wildcards.sample} > {output}
#         '''


rule cutadapt_parse_lens:
    input:
#        samplesdir + "/{sample}/log/cutadapt.log",
        lambda wildcards: get_trim_log(LIBTYPES, wildcards.sample),
    output:
        samplesdir + "/{sample}/log/cutadapt.len.tsv",
    script:
        "../scripts/parse-cutadapt-lens.py"


# rule cutadapt_plot_lens:
#     input:
#         samplesdir + "/{sample}/log/cutadapt.len.tsv",
#     output:
#         samplesdir + "/{sample}/log/cutadapt.len.pdf",
#     shell:
#         '''
#         {activ_conda}

#         workflow/scripts/plot-cutadapt-lens.R {input} {output}
#         '''


rule cutadapt_plot_lens:
    input:
        samplesdir + "/{sample}/log/cutadapt.len.tsv",
    output:
        samplesdir + "/{sample}/log/cutadapt.len.pdf",
    script:
        "../scripts/plot-cutadapt-lens.R"


### Mapping statistics
# rule mapped_transcriptome:
#     input:
#         samplesdir + "/{sample}/db/sqlite.db",
# #        samplesdir + "/{sample}/db/annot_trans_5tera.done",
#         samplesdir + "/{sample}/db/annot_trans_adapt.done",
#         samplesdir + "/{sample}/db/annot_genome_adapt.done",
#     output:
#         samplesdir + "/{sample}/log/sqldb.mapping-stats.transcriptome.txt",
#         samplesdir + "/{sample}/db/sqldb.mapping-stats.transcriptome.done",
#     params:
#         libtype=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample),
#     shell:
#         '''
#         {activ_conda}

#         workflow/scripts/mapping-stats-transcriptome.sh {input[0]} {params.libtype} > {output[0]} && \
#         touch {output[1]}
#         '''


rule mapped_transcriptome:
    input:
        samplesdir + "/{sample}/db/sqlite.db",
#        samplesdir + "/{sample}/db/annot_trans_5tera.done",
        # samplesdir + "/{sample}/db/annot_trans_adapt.done",
        # samplesdir + "/{sample}/db/annot_genome_adapt.done",
        samplesdir + "/{sample}/db/annot_sqldb.done"
    output:
        stats = samplesdir + "/{sample}/log/sqldb.mapping-stats.transcriptome.txt",
        done = samplesdir + "/{sample}/db/sqldb.mapping-stats.transcriptome.done",
    params:
        libtype=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample),
#    script:
#        "../scripts/mapping-stats-transcriptome.sh"
    shell:
        '''
        {activ_conda}

        adapt_cols={params.libtype}

        # For 5tera3, iterate over both adapters
        if [ {params.libtype} == "5tera3" ]; then
            adapt_cols="rel5 rel3"
        elif [ {params.libtype} == "tera3" ]; then
            adapt_cols="rel3"
        elif [ {params.libtype} == "5tera" ]; then
            adapt_cols="rel5"
        fi

        echo "All reads mapped to protein-coding transcripts (primary)" > {output.stats}
        sqlite3 {input[0]} 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> {output.stats}
        echo "Unambiguously mapped reads to protein-coding transcripts" >> {output.stats}
        sqlite3 {input[0]} 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND (best_hit == 1 AND number_of_best_hits == 1));' >> {output.stats}

        for adapt_col in ${{adapt_cols}}; do
            echo "Mapped reads to protein-coding transcripts without ${{adapt_col}} adapter (primary)" >> {output.stats}
            sqlite3 {input[0]} "SELECT COUNT (DISTINCT qname) FROM transcr WHERE (${{adapt_col}} IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));" >> {output.stats}
            echo "Mapped reads to protein-coding transcripts with ${{adapt_col}} adapter (primary)" >> {output.stats}
            sqlite3 {input[0]} "SELECT COUNT (DISTINCT qname) FROM transcr WHERE (${{adapt_col}} IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));" >> {output.stats}

            echo "All mapped protein-coding transcripts without ${{adapt_col}} adapter (from primary mappings)" >> {output.stats}
            sqlite3 {input[0]} "SELECT COUNT (DISTINCT rname) FROM transcr WHERE (${{adapt_col}} IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));" >> {output.stats}
            echo "All mapped protein-coding transcripts with ${{adapt_col}} adapter (from primary mappings)" >> {output.stats}
            sqlite3 {input[0]} "SELECT COUNT (DISTINCT rname) FROM transcr WHERE (${{adapt_col}} IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 16) == 0) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));" >> {output.stats}
        done

        if [ `wc -l {output.stats} | cut -d ' ' -f1` -eq 12 ] || [ `wc -l {output.stats} | cut -d ' ' -f1` -eq 20 ] ; then
            touch {output.done}
        fi
        '''

# rule mapped_genome:
#     input:
#         samplesdir + "/{sample}/db/sqlite.db",
#         samplesdir + "/{sample}/db/sqldb.mapping-stats.transcriptome.done",
# #        samplesdir + "/{sample}/db/annot_genome_5tera.done",
#         samplesdir + "/{sample}/db/annot_trans_adapt.done",
#         samplesdir + "/{sample}/db/annot_genome_adapt.done",
#     output:
#         samplesdir + "/{sample}/log/sqldb.mapping-stats.genome.txt",
#         samplesdir + "/{sample}/db/sqldb.mapping-stats.genome.done",
#     shell:
#         '''
#         {activ_conda}

#         workflow/scripts/mapping-stats-genome.sh {input[0]} > {output[0]} && \
#         touch {output[1]}
#         '''

rule mapped_genome:
    input:
        samplesdir + "/{sample}/db/sqlite.db",
        samplesdir + "/{sample}/db/sqldb.mapping-stats.transcriptome.done",
        # samplesdir + "/{sample}/db/annot_trans_adapt.done",
        # samplesdir + "/{sample}/db/annot_genome_adapt.done",
        samplesdir + "/{sample}/db/annot_sqldb.done"
    output:
        stats=samplesdir + "/{sample}/log/sqldb.mapping-stats.genome.txt",
        done=samplesdir + "/{sample}/db/sqldb.mapping-stats.genome.done",
#    script:
#        "../scripts/mapping-stats-genome.sh"
    shell:
        '''
        {activ_conda}

        echo "All reads mapped to genome (primary)" > {output.stats}
        sqlite3 {input[0]} 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0));' >> {output.stats}
        echo "Unambiguously mapped reads to genome (unambiguous, primary flag not considered)" >> {output.stats}
        sqlite3 {input[0]} 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (best_hit == 1 AND number_of_best_hits == 1));' >> {output.stats}

        if [ `wc -l {output.stats} | cut -d ' ' -f1` -eq 4 ]; then
            touch {output.done}
        fi
        '''


### Read lengths
# Trimmed fastq lengths - histogram (read-len-hist.tsv)
rule read_length_hist_all:
    input:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.txt",
    params:
        mem_gb=3,
        tmpdir=lambda wildcards, input: Path(input[0]).parent
    threads: 16
    shell:
        '''
        {activ_conda}

        bbmap=`find ${{CONDA_PREFIX}}/opt/bbmap* -type d -name 'current'`

        java -ea -Djava.io.tmpdir={params.tmpdir} -Xmx{params.mem_gb}g -cp $bbmap jgi.MakeLengthHistogram \
            in={input} bin=1 nzo=f qin=33 out={output}
        '''


# rule read_length_hist_all_parse:
#     input:
#         samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.txt",
#     output:
#         samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.tsv",
#     params:
#         keepzero=False,
#     shell:
#         '''
#         python3 workflow/scripts/parse-bbmap-lens.py {input} {wildcards.sample} {params.keepzero} \
#             > {output}
#         '''


rule read_length_hist_all_parse:
    input:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.tsv",
    params:
        keepzero=False,
    script:
        "../scripts/parse-bbmap-lens.py"


# Trimmed fastq lengths - per read (read-len.tsv)
# rule read_length_all:
#     input:
#         samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
#     output:
#         samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz",
#     shell:
#         '''
#         {activ_conda}
#
#         python3 workflow/scripts/fastq-read-length.py {input} | gzip -c > {output}
#         '''


rule read_length_all:
    input:
        samplesdir + "/{sample}/fastq/reads.1.sanitize.adapt_trim.fastq.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz",
    script:
        "../scripts/fastq-read-length.py"


# Mapped total length (read-len.tsv) and only the mapped portion (mapped-len.tsv)
rule mapped_length_genome:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.tsv.gz",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz",
    params:
        flags="-F 4 -F 256 -F 2048",
        tmpdir=lambda wildcards, input: Path(input[0]).parent
    threads: 16
    shell:
        '''
        {activ_conda}

        # Aligned read length - includes softclipping/extensive splicing
        echo -e "read_id\\tchrom\\tlength" | gzip -c > {output.read_len} && \
        samtools view -@ {threads} {params.flags} {input} \
            | awk 'BEGIN{{OFS="\\t"}}{{print $1,$3,length($10)}}' \
            | sort -T {params.tmpdir} --parallel={threads} -k3,3nr \
            | gzip -c >> {output.read_len}

        # Aligned section of the read (sum of bases); counts M, X, = and D, doesn't count N, I, ... - excludes softclipped bases
        echo -e "read_id\\tchrom\\tlength" | gzip -c > {output.mapped_len} && \
        samtools view -@ {threads} {params.flags} {input} \
            | perl -slane '$l = 0; $F[5] =~ s/(\d+)[MX=D]/$l+=$1/eg; print $F[0],"\\t",$F[2],"\\t",$l' \
            | sort -T {params.tmpdir} --parallel={threads} -k3,3nr \
            | gzip -c >> {output.mapped_len}
        '''


use rule mapped_length_genome as mapped_length_transcriptome with:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.tsv.gz",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.align-len.full.tsv.gz",
    params:
        flags="-F 4 -F 256 -F 2048 -F 16",
        tmpdir=lambda wildcards, input: Path(input[0]).parent,


# Plots
# rule total_aligned_genome_length_plot:
#     input:
#         total_len=samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz",
#         aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz",
#     output:
#         samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.trim-vs-aligned.pdf",
#     params:
#         len_cap=5000,
#     shell:
#         '''
#         {activ_conda}

#         workflow/scripts/plot-aligned-vs-fastq-length.R \
#             {input.total_len} {input.aligned_len} {params.len_cap} \
#             {output}
#         '''


rule total_aligned_genome_length_plot:
    input:
        total_len=samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz",
        aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.trim-vs-aligned.pdf",
    params:
        len_cap=5000,
    script:
        "../scripts/plot-aligned-vs-fastq-length.R"


use rule total_aligned_genome_length_plot as total_aligned_transcriptome_length_plot with:
    input:
        total_len=samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz",
        aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.align-len.full.tsv.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.len.trim-vs-aligned.pdf",


# rule mapped_length_genome_plot:
#     input:
#         mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.tsv.gz",
#         aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz",
#     output:
#         samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.mapped-vs-aligned.pdf",
#     params:
#         len_cap=5000,
#     shell:
#         '''
#         {activ_conda}

#         workflow/scripts/plot-aligned-length.R \
#             {input.mapped_len} {input.aligned_len} {params.len_cap} \
#             {output}
#         '''


rule mapped_length_genome_plot:
    input:
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.tsv.gz",
        aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.mapped-vs-aligned.pdf",
    params:
        len_cap=5000,
    script:
        "../scripts/plot-aligned-length.R"


use rule mapped_length_genome_plot as mapped_length_transcriptome_plot with:
    input:
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.tsv.gz",
        aligned_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.align-len.full.tsv.gz",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.len.mapped-vs-aligned.pdf",


rule mapping_stats_genome:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapping-stats.tsv",        
    params:
        samflags="-F 4 -F 256 -F 2048",
    threads: 16
    shell:
        '''
        echo -e "stat\t{wildcards.sample}" > {output}

        samtools stats -@ {threads} {params.samflags} {input} \
            | grep '^SN' \
            | cut -f 2- \
            | cut -f1,2 \
            >> {output}
        '''


use rule mapping_stats_genome as mapping_stats_transcriptome with:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapping-stats.tsv",        
    params:
        samflags="-F 4 -F 256 -F 2048 -F 16",    