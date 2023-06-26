"""
Add allele length delta columns onto a bed file
"""
import sys
import pysam
import truvari
import numpy as np


vcf = pysam.VariantFile(sys.argv[1])
regions = truvari.opt_gz_open(sys.argv[2])

for reg in regions:
    reg = reg.strip()
    chrom, start, end = reg.split("\t")[:3]
    start = int(start)
    end = int(end)
    d1 = 0 
    d2 = 0
    for entry in vcf.fetch(chrom, start, end):
        st, ed = truvari.entry_boundaries(entry)
        if not (start <= st <= ed <= end):
            continue
        sz = len(entry.ref) - len(entry.alts[0])
        d1 += sz if entry.samples[0]['GT'][0] == 1 else 0
        d2 += sz if entry.samples[0]['GT'][1] == 1 else 0
    print(reg, d1, d2, sep='\t')
