#!/usr/bin/env python3
import os
import sys

#ifile="/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/cutadapt.log"

def cutadapt_parse_lens(ifile, library):
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
    print("library" + "\t" + adapter_lines[0])
    for line in adapter_lines:
        if line[:1].isdigit():
            length = line.split('\t')[0]
            count = line.split('\t')[1]
            expect = line.split('\t')[2]
            max_err = line.split('\t')[3]
            err_counts = line.split('\t')[4]
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(library, length, count, expect, max_err, err_counts))

usage = "Usage: python3 " + os.path.basename(__file__) + " cutadapt.log samplename"

if len(sys.argv) != 3:
    print(usage)
else:
    ifile=sys.argv[1]
    if sys.argv[1] in ("-h --help"):
        print(usage)
    else:
        library=sys.argv[2]
        cutadapt_parse_lens(ifile, library)
