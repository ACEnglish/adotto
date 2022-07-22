#!/usr/bin/env python3
"""
Make simple stats on a bed file - works inside a pipe so reads from stdin and writes the same to stdout
make a file bed_stats.json
"""
import sys
import json
import pandas as pd

suffix_name = '_' + sys.argv[1] if len(sys.argv) > 1 else ""
spans = []
for line in sys.stdin:
    chrom, start, end = line.strip().split('\t')[:3]
    start = int(start)
    end = int(end)
    spans.append(abs(end - start))
    sys.stdout.write(line)

spans = pd.Series(spans, name="span")
desc = spans.describe().astype(int).to_string()
sys.stderr.write(desc + '\n')
with open(f"input_spans{suffix_name}.txt", 'w') as fout:
    fout.write(desc + '\n')
