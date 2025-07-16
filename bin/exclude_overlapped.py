#!/usr/bin/env python3

import sys

lengths_file = sys.argv[1]
contigs_file = sys.argv[2]
exclude_file = lengths_file.replace('_lengths.tsv', '_exclude.txt')
with open(lengths_file, 'r') as f1, open(contigs_file, 'r') as f2, open(exclude_file, 'w') as f3:
    lengths = {}
    for line in f2:
        contig_id, length = line.split()[0:2]
        lengths[contig_id] = int(length)
    for line in f1:
        contig_id, length = line.split()[0:2]
        length = int(length)
        coverage = length / lengths[contig_id] * 100
        if coverage > 90:  # Adjust the threshold as needed
            f3.write(contig_id + '\n')