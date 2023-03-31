from pathlib import Path

#def get_links(assembly):
    # if assembly == "hg38":
    #     return {'org':"Homo sapiens",
    #             'link':"ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz",
    #             'link_gtf':"ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz",
    #             'gtrna':"http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz",
    #             'gtrna_bed':datadir + "/{assembly}/hg38-tRNAs.bed"}
    # elif assembly == "mm10":
    #     return {'org':"Mus musculus",
    #             'link':"ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",
    #             'link_gft':"ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz"}

### SILVA rRNA
rule silva_download:
    output:
#        silva_lsu=temp(datadir + "/silva/SILVA_132_LSURef_tax_silva_trunc.fasta.gz"),
#        silva_ssu=temp(datadir + "/silva/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz")
       silva_lsu=datadir + "/silva/SILVA_132_LSURef_tax_silva_trunc.fasta.gz",
       silva_ssu=datadir + "/silva/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz",
    params:
#        link_lsu="https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz",
#        link_ssu="https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz",
        link_lsu=REF_LINKS['silva']['lsu'],
        link_ssu=REF_LINKS["silva"]['ssu'],
    shell:
        '''
        wget \
             {params.link_lsu} \
            -O {output.silva_lsu}
        wget \
             {params.link_ssu}\
            -O {output.silva_ssu}
        '''

rule silva_merge:
    input:
        silva_lsu=datadir + "/silva/SILVA_132_LSURef_tax_silva_trunc.fasta.gz",
        silva_ssu=datadir + "/silva/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz",
    output:
        datadir + "/silva/ribosomal.fa"
    shell:
        '''
        {activ_perl}

        # Merge all ribosomal sequences together eliminating duplicates
        zcat \
            {input.silva_lsu} \
            {input.silva_ssu} \
            | fasta-rm-dup-seqs \
            | fasta-unique-names \
            | sed '/^[^>]/ y/uU/tT/' \
            > {output}
        '''

rule silva_sanitize:
    input:
        datadir + "/silva/ribosomal.fa"
    output:
        datadir + "/silva/ribosomal.bed"
    shell:
        '''
        {activ_perl}

        # Create simple BED file where stop position corresponds to sequence length.
        fasta-sanitize-header \
            --input {input} \
            --delim : \
            --keep 0 \
            | fasta-to-sizes-bed \
                --name-suffix ribo \
                > {output}
        '''

rule silva_extract:
    input:
        datadir + "/silva/ribosomal.fa"
    output:
        datadir + "/{assembly}/ribosomal.fa"
    params:
#        org="Homo sapiens"
#        org=lambda wildcards: get_links(wildcards.assembly)['org']
        org=lambda wildcards: REF_LINKS[wildcards.assembly]['org']
    shell:
        '''
        {activ_perl}

        # Isolate organism ribosomal genes and eliminate superfluous info from header.
        grep -A 1 --no-group-separator "{params.org}" {input} \
            | fasta-sanitize-header \
                --input - \
                --delim : \
                --keep 0 \
                > {output}
        '''

### GtRNAdb tRNA
rule gtrna_download:
    output:
        tar=datadir + "/{assembly}/tRNAs.tar.gz",
        bed=datadir + "/{assembly}/{assembly}-tRNAs.bed"
    params:
#        link=lambda wildcards: get_links(wildcards.assembly)['gtrna']
        link=lambda wildcards: REF_LINKS[wildcards.assembly]['gtrna'],
        datadir=datadir,
    shell:
        '''
        wget {params.link} -O {output.tar}
        tar -xvzf {output.tar} -C "{params.datadir}/{wildcards.assembly}"
        '''

# TODO: Make the chrosome conversion universal; right now (03/22/2023) it is specific for human and it should work for mouse
#           Use https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2ensembl.txt -O UCSC2ensembl.txt conversions and then replace column Python script (have to find that one)
# TODO: Decide whether use resources: tmpdir instead of params: tmpdir
rule gtrna_extract:
    input:
        datadir + "/{assembly}/{assembly}-tRNAs.bed",
    output:
        datadir + "/{assembly}/tRNA.bed",
    params:
        tmpdir=lambda wildcards, input: Path(input[0]).parent
    threads: 4
    shell:
        '''
        cat {input} \
            | sed 's/^chrM/MT/g' \
            | sed 's/^chr//g' \
            | sed 's/chr1_KI270713v1_random/KI270713.1/g' \
            | sort --parallel={threads} -T {params.tmpdir} -k1,1 -k2,2 | cut -f1-6 > {output}
        '''

rule trna_rrna_merge:
    input:
        trna=datadir + "/{assembly}/tRNA.bed",
        rrna=datadir + "/{assembly}/rRNA.bed",
    output:
        datadir + "/{assembly}/rRNA_tRNA.bed"
    shell:
        '''
        cat {input.rrna} {input.trna} > {output}
        '''

### Genome
rule genome_download:
    output:
        datadir + "/{assembly}/genome.fa",
    params:
#        link="ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz",
#        link=lambda wildcards: get_links(wildcards.assembly)['genome'],
        link=lambda wildcards: REF_LINKS[wildcards.assembly]['genome'],
    shell:
        '''
        {activ_perl}

        wget -qO- {params.link} \
            | gunzip -c \
            | clean-genome-headers --fasta - \
            > {output}
        '''

rule genome_index:
    input:
        datadir + "/{assembly}/genome.fa"
    output:
        fai=datadir + "/{assembly}/genome.fa.fai",
        sizes=datadir + "/{assembly}/chrom.sizes",
    shell:
        '''
        samtools faidx {input}
        cut -f1-2 {output.fai} > {output.sizes}
        '''

### Gene annotation
rule annotation_download:
    output:
        temp(datadir + "/{assembly}/ensembl_genes.orig.gtf")
    params:
#        link_gtf="ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz"
#        link=lambda wildcards: get_links(wildcards.assembly)['gtf']
        link=lambda wildcards: REF_LINKS[wildcards.assembly]['gtf'],
        datadir=datadir,
    shell:
        '''
        wget -qO- {params.link} \
            | gunzip -c \
            > {output}
        ln -sf $(basename {output}) {params.datadir}/{wildcards.assembly}/$(basename {params.link} .gz)
        '''

rule annotation_clean:
    input:
        datadir + "/{assembly}/ensembl_genes.orig.gtf"
    output:
        datadir + "/{assembly}/genes-polya.gtf",
        datadir + "/{assembly}/genes.gtf",
        datadir + "/{assembly}/genes-total.gtf",
    shell:
        '''
        {activ_perl}

        cat {input} \
            | clean-gtf-lines-polya --gtf - \
            > {output[0]}
        ln -sf $(basename {output[0]}) {output[1]}

        cat {input} \
            | clean-gtf-lines-total --gtf - \
            > {output[2]}
        '''

rule annotation_elements:
    input:
        datadir + "/{assembly}/genes.gtf"
    output:
        datadir + "/{assembly}/genic_elements.bed",
        datadir + "/{assembly}/genic_elements-total.bed",
    shell:
        '''
        gff-to-genic-elements-bed \
            --input {input} \
            > {output[0]}
        ln -sf $(basename {output[0]}) {output[1]}
        '''

rule annotation_elements_extract:
    input:
        datadir + "/{assembly}/genic_elements.bed",
        datadir + "/{assembly}/genic_elements-total.bed",
    output:
        multiext(datadir + "/{assembly}/genic_elements", ".utr5.bed", ".utr3.bed", ".cds.bed", ".ncrna.bed", ".mrna.bed"),
        multiext(datadir + "/{assembly}/genic_elements-total", ".utr5.bed", ".utr3.bed", ".cds.bed", ".ncrna.bed", ".mrna.bed"),
    params:
        datadir=datadir,
    shell:
        '''
        # Create separate files for each genic element
        for element in "utr5" "utr3" "cds" "ncrna" "mrna"; do
            if grep -P -q ":${{element}}\t" {input}; then
                grep -P ":${{element}}\t" {input} > {params.datadir}/{wildcards.assembly}/genic_elements.${{element}}.bed
            else
                touch {params.datadir}/{wildcards.assembly}/genic_elements.${{element}}.bed
            fi

            # Create separate files for each genic element - total RNA - just copy the same thing as for polyA
            ln -sf genic_elements.${{element}}.bed {params.datadir}/{wildcards.assembly}/genic_elements-total.${{element}}.bed
        done
        '''

rule transcripts_extract:
    input:
        genome=datadir + "/{assembly}/genome.fa",
        gtf=datadir + "/{assembly}/genes.gtf",
        gtf_total=datadir + "/{assembly}/genes-total.gtf",
        gtf_ens=datadir + "/{assembly}/ensembl_genes.orig.gtf",
    output:
        trans=datadir + "/{assembly}/transcripts.fa",
        trans_total=datadir + "/{assembly}/transcripts-total.fa",
        gtf=datadir + "/{assembly}/ensembl_transcripts_protein_coding.gtf",
    shell:
        '''
        {activ_conda}

        gffread -w {output.trans} -g {input.genome} {input.gtf} # Poly(A)
        gffread -w {output.trans_total} -g {input.genome} {input.gtf_total} # Total

        cat {input.gtf_ens} | grep "transcript_biotype \\"protein_coding\\"" > {output.gtf}
        '''

# TODO: Change {output.index) to point to all the gmap index file and to the directory
rule gtf_map_rrna:
    input:
        genome=datadir + "/{assembly}/genome.fa",
        rrna=datadir + "/{assembly}/ribosomal.fa",
    output:
        index=directory(datadir + "/{assembly}/gmap-2019-09-12"),
        gff=datadir + "/{assembly}/ribosomal.gmap.gff3",
        summary=datadir + "/{assembly}/ribosomal.gmap-summary.txt",
        gtf=datadir + "/{assembly}/ribosomal.gmap.gtf",
    threads: 16
    shell:
        '''
        {activ_conda}

        mkdir -p {output.index}

        gmap_build -d genome -D {output.index} {input.genome} # Make index

        gmap -d genome -D {output.index} {input.rrna} -t {threads} \
            --suboptimal-score=0.0 -f gff3_match_cdna \
            | sed "s/\tcDNA_match\t/\texon\t/g" | sed "s/\tgenome\t/\tsilva\t/g" \
            > {output.gff} # Map ribosomal RNA to genome and get gff3 output to be added to the gene annotation
        gmap -d genome -D {output.index} {input.rrna} -t {threads} \
            --suboptimal-score=0.0 -S > {output.summary} # Map ribosomal RNA to genome and get txt output for manual check

        gff2gtf-gmap {output.gff} \
            | sed "s/\tgmapidx\t/\tsilva\t/g" \
            | sed "s/;$/; gene_biotype \"rRNA\"; transcript_biotype \"rRNA\";/g" \
            | sort -k1,1 -k4,4n > {output.gtf} # Convert GMAP gff3 to gtf
        '''

rule gtf_add_rrna:
    input:
        gtf_rrna=datadir + "/{assembly}/ribosomal.gmap.gtf",
        gtf=datadir + "/{assembly}/ensembl_genes.orig.gtf",
        genome=datadir + "/{assembly}/genome.fa",
    output:
        bed=datadir + "/{assembly}/rRNA.bed",
        trans_wrrna=datadir + "/{assembly}/ensembl-transcripts-wRibo.fa",
        gtf_worrna=temp(datadir + "/{assembly}/ensembl_genes.gtf.tmp"),
        gtf_wrrna=temp(datadir + "/{assembly}/ensembl_genes.gtf.tmp2"),
        gtf=datadir + "/{assembly}/ensembl_genes.gtf",
        gtf_gz=datadir + "/{assembly}/ensembl_genes.gtf.gz",
    params:
#        link="ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz"
#        link=lambda wildcards: get_links(wildcards.assembly)['gtf']
        link=lambda wildcards: REF_LINKS[wildcards.assembly]['gtf']
    shell:
        '''
        {activ_conda}

        gtf2bed6 {input.gtf_rrna} | cut -f1-6 > {output.bed}

        cat {input.gtf} | grep -v " \"rRNA\";" > {output.gtf_worrna} # Remove annotated rRNA from Ensembl but keep ribosomal - SILVA rRNA db doesn't annotate those
        cat {output.gtf_worrna} {input.gtf_rrna} > {output.gtf_wrrna}
        (grep "^#" {output.gtf_wrrna}; grep -v "^#" {output.gtf_wrrna} | sort -k1,1 -k4,4n) > {output.gtf} # Add SILVA rRNA to Ensembl
        gzip -c {output.gtf} > {output.gtf_gz}
        ln -sf {output.gtf_gz} $(dirname {output.gtf})/$(basename {params.link})

        gffread -w {output.trans_wrrna} -g {input.genome} {output.gtf} # All the transcripts with rRNA
        '''

### Minimap2
rule index_genome_minimap2:
    input:
        datadir + "/{assembly}/genome.fa",
    output:
        datadir + "/{assembly}/minimap2.17/genome.k12.mmi",
    params:
        k=12,
        odir=lambda wildcards, output: Path(output[0]).parent,
    shell:
        '''
        {activ_conda}

        if [ ! -d "{params.odir}" ]; then
            mkdir -p {params.odir}
        fi

        minimap2 -k {params.k} -d {output} {input}
        '''

use rule index_genome_minimap2 as index_trans_polya_minimap2 with:
    input:
        trans=datadir + "/{assembly}/transcripts.fa",
    output:
        mmi_trans=datadir + "/{assembly}/minimap2.17/transcripts.k12.mmi",

use rule index_genome_minimap2 as index_trans_total_minimap2 with:
    input:
        trans_total=datadir + "/{assembly}/transcripts-total.fa",
    output:
        mmi_total=datadir + "/{assembly}/minimap2.17/transcripts-total.k12.mmi",

use rule index_genome_minimap2 as index_trans_wrrna_minimap2 with:
    input:
        trans_wrrna=datadir + "/{assembly}/ensembl-transcripts-wRibo.fa",
    output:
        mmi_wribo=datadir + "/{assembly}/minimap2.17/ensembl-transcripts-wRibo.k12.mmi",


### STAR
# Note: For human-sized genomes, STAR will need ~33 GB RAM unless you change --genomeSAsparseD to 2 (default: 1) which will fit the genome to ~16 GM RAM but will make the search slower   input:
rule index_woannot_star:
    input:
        genome=datadir + "/{assembly}/genome.fa",
    output:
        index=datadir + "/{assembly}/STAR-2.7.2b/SAindex",
    params:
#        tmpdir=str(lambda wildcards, output: Path(output[0]).parent) + "_STARtmp_woannot", # This won't work because lambda is evaluated in the shell: part (?)
        odir=lambda wildcards, output: Path(output[0]),
        genomeSAsparseD=1,
        mem_b=31000000000,
    threads: 32
    shell:
        '''
        {activ_conda}

        mkdir -p {params.odir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeSAsparseD {params.genomeSAsparseD} \
            --limitGenomeGenerateRAM {params.mem_b} \
            --outTmpDir {params.odir}/_STARtmp \
            --outFileNamePrefix {params.odir} \
            --genomeDir {params.odir}/ \
            --genomeFastaFiles {input.genome}
        '''

rule index_wannot_star:
    input:
        genome=datadir + "/{assembly}/genome.fa",
        gtf=datadir + "/{assembly}/ensembl_genes.gtf",
    output:
        index=datadir + "/{assembly}/STAR-2.7.2b-annot/SAindex",
    params:
#        tmpdir=str(lambda wildcards, output: Path(output[0]).parent) + "_STARtmp_wannot", # This won't work because lambda is evaluated in the shell: part (?)
        odir=lambda wildcards, output: Path(output[0]),
        genomeSAsparseD=1,
        mem_b=31000000000,
    threads: 32
    shell:
        '''
        {activ_conda}

        mkdir -p {params.odir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeSAsparseD {params.genomeSAsparseD} \
            --limitGenomeGenerateRAM {params.mem_b} \
            --outTmpDir {params.odir}/_STARtmp \
            --outFileNamePrefix {params.odir} \
            --genomeDir {params.odir}/ \
            --genomeFastaFiles  {input.genome} \
            --sjdbGTFfile  {input.gtf} \
            --sjdbOverhang 100
        '''
