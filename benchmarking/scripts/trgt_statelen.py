"""
Parse a bench directory and report the len(REF) len(ALT) and state of files
"""
import os
import sys
import pysam
import pandas as pd

in_dir = sys.argv[1]

files = [('fn', 'fn.vcf.gz'),
         ('fp', 'fp.vcf.gz'),
         ('tp', 'tp-call.vcf.gz')]

rows = []
for state, i in files:
    v = pysam.VariantFile(os.path.join(in_dir, i))
    for entry in v:
        rows.append([state, len(entry.ref), len(entry.alts[0])])
data = pd.DataFrame(rows, columns=["state", "ref", "alt"])
data.to_csv('/dev/stdout', sep="\t", index=False)
