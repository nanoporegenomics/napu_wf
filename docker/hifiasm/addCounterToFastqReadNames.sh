#!/usr/bin/env bash

# addCounterToFastqReadNames.sh
#
# Fixes duplicate read names in a FASTQ file by appending a global counter
# to every read name, guaranteeing uniqueness.
#
# Usage:
#   ./addCounterToFastqReadNames.sh input.fastq > output.fastq
#   ./addCounterToFastqReadNames.sh input.fastq.gz > output.fastq
#   ./addCounterToFastqReadNames.sh input.fastq.gz | gzip > output.fastq.gz


set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $(basename "$0") <input.fastq[.gz]>" >&2
  exit 1
fi

INPUT="$1"

# plain or gzipped input 
if [[ "$INPUT" == *.gz ]]; then
  READER="gzip -dc"
else
  READER="cat"
fi

# Rewrite read names 
# A FASTQ record is exactly 4 lines:
#   1. @<read_name> 
#   2. sequence
#   3. + 
#   4. quality scores
#
# For each record:
#   - Line 1: strip everything after the first space, append _<counter>; Dont 
#               worry about methylation tags being carried over, this is just for 
#               downsampling prior to assembly. BAM --> fastQ --> renamedIDs --> hifiasm
#   - Line 3: normalise to bare "+" 
#   - Lines 2 & 4: printed unchanged

$READER "$INPUT" | awk '
BEGIN { counter = 0 }

# Line 1 of record: the @ header
(NR % 4 == 1) {
    counter++
    # Grab just the read name (up to first space), drop any existing comment
    split($0, parts, " ")
    name = parts[1]          
    # Strip leading @ for manipulation, then reattach
    sub(/^@/, "", name)
    print "@" name "_" counter
    next
}

# Line 3 of record: the + separator 
(NR % 4 == 3) {
    print "+"
    next
}

# Lines 2 and 4: sequence and quality — print unchanged
{ print }
'
