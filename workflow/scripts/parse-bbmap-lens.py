#!/usr/bin/env python3
import os
import sys

#filename="/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.adapt_trim.read-len.tsv"

def bbmap_parse_lens(filename, library, keepzero=False):
    # Read filename line by lin
    with open(filename, 'r') if filename != "-" else sys.stdin as f:
        content = f.read()

    # Find the adapter trimming section in the output
    start_index = content.find('#Length\t')
    end_index = len(content)

    # Extract the adapter trimming information
    len_info = content[start_index:end_index]

    # Split the adapter information into lines
    len_lines = len_info.split('\n')

    # Extract the length and number of trimmed adapters for each adapter sequence
    print('{}\t{}\t{}'.format("library", "length", "count")) #  len_lines[0])
    for line in len_lines:
        if line[:1].isdigit():
            length = line.split('\t')[0]
            count = line.split('\t')[1]            

            if int(count) > 0 or keepzero == "True":
                print('{}\t{}\t{}'.format(library, length, count))

usage = "Usage: python3 " + os.path.basename(__file__) + " bbmap.read-len.txt samplename keepzero[True|False; default: False]. Use \"-\" for stdin."

if len(sys.argv) < 3:
    print(usage)
else:
    filename=sys.argv[1]
#    if sys.argv[1] in ("-h --help"):
#        print(usage)
#    else:
    library=sys.argv[2]
    
    keepzero=False
    if len(sys.argv) == 4:
        keepzero=sys.argv[3]

    bbmap_parse_lens(filename, library, keepzero)
