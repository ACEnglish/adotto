import sys
import pysam
import joblib
import truvari
import numpy as np

bed_fn = sys.argv[1]
vcf_fn = sys.argv[2]
BUFFERS = range(0, 201, 10)

bed = truvari.opt_gz_open(bed_fn)
vcf = pysam.VariantFile(vcf_fn)

wp = np.zeros((len(BUFFERS), 2))
for line in bed:
    data = line.strip().split('\t')
    chrom, start, end = data[:3]
    start = int(start)
    end = int(end)
    for entry in vcf.fetch(chrom, start, end):
        sz = truvari.entry_size(entry)
        if sz < 5:
            continue
        s,e = truvari.entry_boundaries(entry)
        for p, b in enumerate(BUFFERS):
            if start - b <= s <= e <= end + b:
                wp[p][0] += 1
            elif start - b<= s <= end or start <= e <= end + b:
                wp[p][1] += 1

joblib.dump(wp, 'wp.jl')

