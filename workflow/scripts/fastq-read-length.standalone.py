#!/usr/bin/env python3
#
# Get tab separated read length per read from fastq.gz/fastq 
#
# One argument - input fastq.gz or fastq if stdin 
#
import os
import sys
import gzip

def parse_fastq(filename):
    with gzip.open(filename, 'rb') if filename != "-" else sys.stdin as f:
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


if len(sys.argv) != 2:
    print(usage)
else:
    filename=sys.argv[1]

    #filename = '/home/jan/playground/TERA-Seq_snakemake/data/samples/sc.test.1/fastq/reads.1.fastq.gz'
    print("read_id" + "\t" + "length")
    for name, length in parse_fastq(filename):
        print(f"{name}\t{length}")