import sys
import pysam
import truvari
import random
import math

bed = sys.argv[1]
vcf = sys.argv[2]
sample = sys.argv[3]

# for each region, calculate the maximum allele delta 

vcf = pysam.VariantFile(vcf)
fh = truvari.opt_gz_open(bed)
for line in fh:
    line = line.strip()
    chrom, start, end = line.split('\t')[:3]
    start = int(start)
    end = int(end)
    h1 = 0
    h2 = 0
    for entry in vcf.fetch(chrom, start, end):
        if 1 not in entry.samples[sample]["GT"]:
            continue
        st, ed = truvari.entry_boundaries(entry)
        if not (start <= st and ed <= end):
            continue
        sz = truvari.entry_size(entry)
        ty = truvari.entry_variant_type(entry)
        if ty == truvari.SV.DEL:
            sz = -sz
        if entry.samples[sample]["GT"][0] == 1:
            h1 += sz
        if entry.samples[sample]["GT"][1] == 1:
            h2 += sz
    m_size = None
    if abs(h1) == abs(h2):
        # might have equal length expansion/contraction
        m_size = [h1, h2][round(random.random())]
    if abs(h1) > abs(h2):
        m_size = h1
    if abs(h1) < abs(h2):
        m_size = h2
    print(chrom, start, end, m_size, sep='\t')

