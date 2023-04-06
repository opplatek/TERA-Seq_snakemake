#!/bin/bash
#
# Summarize number of mapped reads to genome reference
#
# Two arguments: sqlite db file, library type [5tera, tera3]
#
# Note: 5tera3 is not implemented as of 04/06/2023
#

db=${snakemake_input[0]}
libtype=${snakemake_params[libtype]}
ofile=${snakemake_output[0]}
done=${snakemake_output[1]}

libtype=`echo $libtype | tr '[:upper:]' '[:lower:]'`

echo "All reads mapped to protein-coding transcripts (primary)" > ${ofile}
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
echo "Unambiguously mapped reads to protein-coding transcripts" >> ${ofile}
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND (best_hit == 1 AND number_of_best_hits == 1));' >> ${ofile}

if [ "$libtype" = "5tera" ]; then
    echo "Mapped reads to protein-coding transcripts without REL5 adapter (primary)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
    echo "Mapped reads to protein-coding transcripts with REL5 adapter (primary)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}

    echo "All mapped protein-coding transcripts without REL5 adapter (from primary mappings)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
    echo "All mapped protein-coding transcripts with REL5 adapter (from primary mappings)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
elif [ "$libtype" = "tera3" ]; then
    echo "Mapped reads to protein-coding transcripts without REL3 adapter (primary)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
    echo "Mapped reads to protein-coding transcripts with REL3 adapter (primary)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}

    echo "All mapped protein-coding transcripts without REL3 adapter (from primary mappings)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
    echo "All mapped protein-coding transcripts with REL3 adapter (from primary mappings)" >> ${ofile}
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND ((flag & 16) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' >> ${ofile}
fi

if [ `wc -l ${ofile} | cut -d ' ' -f1` -eq 12 ]; then
    touch ${done}
fi