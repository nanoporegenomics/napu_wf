#!/usr/bin/env python3

import sys

mode = sys.argv[1]
length_threshold = int(sys.argv[2])

for line in sys.stdin:
    if line.startswith("#"):
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
