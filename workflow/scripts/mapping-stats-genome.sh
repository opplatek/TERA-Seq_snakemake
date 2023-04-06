#!/bin/bash
#
# Summarize number of mapped reads to genome reference
#
# One argument: sqlite db file
#

db=${snakemake_input[0]}
ofile=${snakemake_output[0]}
done=${snakemake_output[1]}

# snakemake_activ_conda {activ_conda}

echo "All reads mapped to genome (primary)" > ${ofile}
sqlite3 ${db} 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0));' >> ${ofile}
echo "Unambiguously mapped reads to genome (unambiguous, primary flag not considered)" >> ${ofile}
sqlite3 ${db} 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (best_hit == 1 AND number_of_best_hits == 1));' >> ${ofile}

if [ `wc -l ${ofile} | cut -d ' ' -f1` -eq 4 ]; then
    touch ${done}
fi