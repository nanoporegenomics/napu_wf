#!/usr/bin/env python3

import sys

mode = sys.argv[1]
length_threshold = int(sys.argv[2])
sample_name = sys.argv[3]

for line in sys.stdin:
    if line.startswith("#"):
        if line.startswith('#CHROM'):
            tokens = line.strip().split('\t')
            if tokens[-1]=='Sample':
                tokens[-1] = sample_name
                sys.stdout.write('\t'.join(tokens))
        else:
            sys.stdout.write(line)
        continue
    fields = line.split()
    sv_len = abs(len(fields[3]) - len(fields[4]))

    if mode == "greater":
        if sv_len >= length_threshold:
            sys.stdout.write(line)

    if mode == "less":
        if sv_len < length_threshold:
            sys.stdout.write(line)
