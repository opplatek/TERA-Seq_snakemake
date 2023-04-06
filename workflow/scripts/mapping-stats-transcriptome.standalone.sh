#!/bin/bash
#
# Summarize number of mapped reads to genome reference
#
# Two arguments: sqlite db file, library type [5tera, tera3, 5tera3]
#

db=$1
libtype=$2

libtype=`echo $libtype | tr '[:upper:]' '[:lower:]'`

echo "All reads mapped to protein-coding transcripts (primary)"
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
echo "Unambiguously mapped reads to protein-coding transcripts"
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND (best_hit == 1 AND number_of_best_hits == 1));'

if [ "$libtype" = "5tera" ]; then
    echo "Mapped reads to protein-coding transcripts without REL5 adapter (primary)"
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
    echo "Mapped reads to protein-coding transcripts with REL5 adapter (primary)"
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'

    echo "All mapped protein-coding transcripts without REL5 adapter (from primary mappings)"
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
    echo "All mapped protein-coding transcripts with REL5 adapter (from primary mappings)"
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
elif [ "$libtype" = "tera3" ]; then
    echo "Mapped reads to protein-coding transcripts without REL3 adapter (primary)"
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
    echo "Mapped reads to protein-coding transcripts with REL3 adapter (primary)"
    sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'

    echo "All mapped protein-coding transcripts without REL3 adapter (from primary mappings)"
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
    echo "All mapped protein-coding transcripts with REL3 adapter (from primary mappings)"
    sqlite3 $db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
fi