#!/usr/bin/env python3
#
# Get tab separated read length per read from fastq.gz/fastq 
#
# One argument - input fastq.gz or fastq if stdin 
#
import os
import sys
import gzip

def parse_fastq(ifile):
    with gzip.open(ifile, 'rb') if ifile != "-" else sys.stdin as f:
        while True:
            # Read four lines at a time to extract one fastq entry
            name = f.readline().strip()
            if not name: break
            sequence = f.readline().strip()
            f.readline()
            quality = f.readline().strip()

            # Yield the read name and length
            yield str(name, 'utf-8')[1:], len(sequence)

usage = "Usage: python3 " + os.path.basename(__file__) + " fastq.gz. Use \"-\" for stdin (requires fastq)."

ifile=snakemake.input[0]
ofile=snakemake.output[0]

with gzip.open(ofile, 'wb') as out:
    out.write(("read_id" + "\t" + "length\n").encode())

    for name, length in parse_fastq(ifile):
        out.write(f"{name}\t{length}\n".encode())