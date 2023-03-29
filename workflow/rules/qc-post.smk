rule cutadapt_parse_lens:
    input:
        "data/samples/{sample}/log/cutadapt.log",
    output:
        "data/samples/{sample}/log/cutadapt.lens.tsv",
    shell:
        '''
        {activ_conda}
        
        python3 workflow/scripts/parse-cutadapt-lens.py {input} {wildcards.sample} > {output}
        '''


rule cutadapt_plot_lens:
    input:
        "data/samples/{sample}/log/cutadapt.lens.tsv",
    output:
        "data/samples/{sample}/log/cutadapt.lens.pdf",
    shell:
        '''
        {activ_conda}

        workflow/scripts/plot-cutadapt-lens.R {input} {output}
        '''


rule mapped_genome:
    input:
        "data/samples/{sample}/db/sqlite.db",    
#        "data/samples/{sample}/log/annot_genome_rel5.done",
        "data/samples/{sample}/log/annot_genome_adapt.done",
    output:
        "data/samples/{sample}/log/mapping-stats.genome.txt",
        temp("data/samples/{sample}/log/mapping-stats.genome.done"),
    shell:
        '''
        {activ_conda}

        workflow/scripts/mapping-stats-genome.sh {input[0]} > {output[0]} \
        && touch {output[1]}
        '''


rule mapped_transcriptome:
    input:
        "data/samples/{sample}/db/sqlite.db",    
        "data/samples/{sample}/log/mapping-stats.genome.done",
#        "data/samples/{sample}/log/annot_trans_rel5.done",
        "data/samples/{sample}/log/annot_trans_adapt.done",
    output:
        "data/samples/{sample}/log/mapping-stats.transcriptome.txt",
        "data/samples/{sample}/log/mapping-stats.transcriptome.done",
    params:
        libtype=lambda wildcards: get_libtype(LIBTYPES, wildcards.sample),
    shell:
        '''
        {activ_conda}

        workflow/scripts/mapping-stats-transcriptome.sh {input[0]} {params.libtype} > {output[0]} \
        && touch {output[1]}
        '''