#!/bin/bash
#
# Summarize number of mapped reads to genome reference
#
# One argument: sqlite db file
#

db=$1

echo "All reads mapped to genome (primary)"
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0));'
echo "Unambiguously mapped reads to genome (unambiguous, primary flag not considered)"
sqlite3 $db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (best_hit == 1 AND number_of_best_hits == 1));'