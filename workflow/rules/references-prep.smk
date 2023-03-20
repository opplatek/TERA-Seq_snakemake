rule silva_download:
    output:
        "test.txt"
#        "data/SILVA_132_LSURef_tax_silva_trunc.fasta.gz",
#        "data/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz"
    shell:
        '''
        . /root/miniconda3/bin/activate
        conda activate teraseq
        samtools --version > {output}
        '''
#        wget \
#            https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz \
#            -O {output[0]}
#        wget \
#            https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz \
#            -O {output[1]}
#        '''
