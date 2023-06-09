#!/usr/bin/env python3
import os
import sys

#filename="/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.adapt_trim.read-len.tsv"

def bbmap_parse_lens(ifile, library, ofile, keepzero=False):
    # Read filename line by lin
    with open(ifile, 'r') if ifile != "-" else sys.stdin as f:
        content = f.read()

    # Find the adapter trimming section in the output
    start_index = content.find('#Length\t')
    end_index = len(content)

    # Extract the adapter trimming information
    len_info = content[start_index:end_index]

    # Split the adapter information into lines
    len_lines = len_info.split('\n')

    # Extract the length and number of trimmed adapters for each adapter sequence
    with open(ofile, 'w') as out:
        out.write('{}\t{}\t{}\n'.format("library", "length", "count")) #  len_lines[0])
        for line in len_lines:
            if line[:1].isdigit():
                length = line.split('\t')[0]
                count = line.split('\t')[1]            

                if int(count) > 0 or keepzero == "True":
                    out.write('{}\t{}\t{}\n'.format(library, length, count))

usage = "Usage: python3 " + os.path.basename(__file__) + " bbmap.read-len.txt samplename keepzero[True|False; default: False]. Use \"-\" for stdin."

ifile=snakemake.input[0]
library=snakemake.wildcards.sample
ofile=snakemake.output[0]
keepzero=snakemake.params.keepzero

bbmap_parse_lens(ifile, library, ofile, keepzero)
