import pandas as pd

"""
with open("grch38.simplerepeat_merge.bed",'r') as fh:
    for line in fh:
        d = line.strip().split('\t')
        start = int(d[1])
        end = int(d[2])
        print(f"{d[0]}\t{start - 10_000}\t{start}")
        print(f"{d[0]}\t{end}\t{end + 10_000}")
"""
import truvari
fh = truvari.opt_gz_open("grch38.gene.raw.txt.gz")
next(fh)
for line in fh:
    data = line.strip().split('\t')
    if data[10] == "":
        continue
    chrom = data[1]
    pos = int(data[3]) if data[2] == '+' else int(data[4])
    other = pos - 1_000 if data[2] == '+' else pos + 1_000
    if other < 0: continue
    print(f"{chrom}\t{min(pos, other)}\t{max(pos, other)}")
