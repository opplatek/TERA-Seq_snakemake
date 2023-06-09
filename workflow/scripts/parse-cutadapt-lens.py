#!/usr/bin/env python3
import os
#import sys

#ifile="/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/cutadapt.log"

def cutadapt_parse_lens(ifile, library, ofile):
    # Read ifile line by lin
    with open(ifile, 'r') as f:
        content = f.read()

    # Find the adapter trimming section in the output
    start_index = content.find('length\t')
    end_index = len(content)

    # Extract the adapter trimming information
    adapter_info = content[start_index:end_index]

    # Split the adapter information into lines
    adapter_lines = adapter_info.split('\n')

    # Extract the length and number of trimmed adapters for each adapter sequence
    with open(ofile, 'w') as out:
        out.write("library" + "\t" + adapter_lines[0] + "\n")
        for line in adapter_lines:
            if line[:1].isdigit():
                length = line.split('\t')[0]
                count = line.split('\t')[1]
                expect = line.split('\t')[2]
                max_err = line.split('\t')[3]
                err_counts = line.split('\t')[4]
                out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(library, length, count, expect, max_err, err_counts))

usage = "Usage: python3 " + os.path.basename(__file__) + " cutadapt.log samplename"

ifile=snakemake.input[0]
library=snakemake.wildcards.sample
ofile=snakemake.output[0]

cutadapt_parse_lens(ifile, library, ofile)
