rule silva_download:
    output:
        "data/SILVA_132_LSURef_tax_silva_trunc.fasta.gz",
        "data/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz"
    shell:
        '''
        wget \
            https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz \
            -O {output[0]}
        wget \
            https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz \
            -O {output[1]}
        '''
