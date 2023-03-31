rule cutadapt_parse_lens:
    input:
        samplesdir + "/{sample}/log/cutadapt.log",
    output:
        samplesdir + "/{sample}/log/cutadapt.lens.tsv",
    shell:
        '''
        {activ_conda}

        python3 workflow/scripts/parse-cutadapt-lens.py {input} {wildcards.sample} > {output}
        '''


rule cutadapt_plot_lens:
    input:
        samplesdir + "/{sample}/log/cutadapt.lens.tsv",
    output:
        samplesdir + "/{sample}/log/cutadapt.lens.pdf",
    shell:
        '''
        {activ_conda}

        workflow/scripts/plot-cutadapt-lens.R {input} {output}
        '''

### Mapping statistics
rule mapped_genome:
    input:
        samplesdir + "/{sample}/db/sqlite.db",
#        samplesdir + "/{sample}/log/annot_genome_rel5.done",
        samplesdir + "/{sample}/log/annot_genome_adapt.done",
    output:
        samplesdir + "/{sample}/log/mapping-stats.genome.txt",
        temp(samplesdir + "/{sample}/log/mapping-stats.genome.done"),
    shell:
        '''
        {activ_conda}

        workflow/scripts/mapping-stats-genome.sh {input[0]} > {output[0]} && \
        touch {output[1]}
        '''


rule mapped_transcriptome:
    input:
        samplesdir + "/{sample}/db/sqlite.db",
        samplesdir + "/{sample}/log/mapping-stats.genome.done",
#        samplesdir + "/{sample}/log/annot_trans_rel5.done",
        samplesdir + "/{sample}/log/annot_trans_adapt.done",
    output:
        samplesdir + "/{sample}/log/mapping-stats.transcriptome.txt",
        samplesdir + "/{sample}/log/mapping-stats.transcriptome.done",
    params:
        libtype=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample),
    shell:
        '''
        {activ_conda}

        workflow/scripts/mapping-stats-transcriptome.sh {input[0]} {params.libtype} > {output[0]} && \
        touch {output[1]}
        '''


### Read lengths
# trimmed fastq lengths (read-len.txt)
rule read_length_all:
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


rule read_length_all_parse:
    input:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len-hist.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.txt",
    params:
        keepzero=False,
    shell:
        '''
        python3 workflow/scripts/parse-bbmap-lens.py {input} {wildcards.sample} {params.keepzero} \
            > {output}
        '''


# Mapped total length (read-len.txt) and only the mapped portion (mapped-len.txt)
rule mapped_length_genome:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.toGenome.sorted.bam",
    output:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.txt",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.txt",
    params:
        flags="-F 4 -F 256 -F 2048",
        tmpdir=lambda wildcards, input: Path(input[0]).parent
    threads: 16
    shell:
        '''
        {activ_conda}

        # Aligned read length - includes softclipping/extensive splicing
        echo -e "read_id\\tchrom\\tlength" > {output.read_len} && \
        samtools view -@ {threads} {params.flags} {input} \
            | awk 'BEGIN{{OFS="\\t"}}{{print $1,$3,length($10)}}' \
            | sort -T {params.tmpdir} --parallel={threads} -k3,3nr \
            >> {output.read_len}

        # Aligned section of the read (sum of bases); counts M, X, = and D, doesn't count N, I, ... - excludes softclipped bases
        echo -e "read_id\\tchrom\\tlength" > {output.mapped_len} && \
        samtools view -@ {threads} {params.flags} {input} \
            | perl -slane '$l = 0; $F[5] =~ s/(\d+)[MX=D]/$l+=$1/eg; print $F[0],"\\t",$F[2],"\\t",$l' \
            | sort -T {params.tmpdir} --parallel={threads} -k3,3nr \
            >> {output.mapped_len}
        '''


use rule mapped_length_genome as mapped_length_transcriptome with:
    input:
        samplesdir + "/{sample}/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam",
    output:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.txt",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.align-len.full.txt",
    params:
        flags="-F 4 -F 256 -F 2048 -f 16",
        tmpdir=lambda wildcards, input: Path(input[0]).parent,


# Plots
rule total_mapped_genome_length_plot:
    input:
        total_len=samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.txt",
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.trim-vs-mapped.pdf",
    shell:
        '''
        {activ_conda}

        ./workflow/scripts/plot-aligned-vs-fastq-length.R {input.total_len} {input.read_len} {output}
        '''

use rule total_mapped_genome_length_plot as total_mapped_transcriptome_length_plot with:
    input:
        total_len=samplesdir + "/{sample}/log/reads.1.sanitize.adapt_trim.read-len.txt",
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.len.trim-vs-mapped.pdf",


rule mapped_length_genome_plot:
    input:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.mapped-len.full.txt",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.align-len.full.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.toGenome.len.pdf",
    shell:
        '''
        {activ_conda}

        ./workflow/scripts/plot-aligned-length.R {input.read_len} {input.mapped_len} {output}
        '''

use rule mapped_length_genome_plot as mapped_length_transcriptome_plot with:
    input:
        read_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.mapped-len.full.txt",
        mapped_len=samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.align-len.full.txt",
    output:
        samplesdir + "/{sample}/log/reads.1.sanitize.noribo.toTranscriptome.len.pdf",
