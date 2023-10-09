import sys
import pysam
import truvari
from collections import Counter, defaultdict
vcf_fn, bed_fn = sys.argv[1:]

vcf = pysam.VariantFile(vcf_fn)

bed = truvari.opt_gz_open(bed_fn)
counts = [0, 0]
for line in bed:
    data = line.strip().split('\t')
    chrom = data[0]
    start = int(data[1])
    end = int(data[2])
    for entry in vcf.fetch(chrom, start, end):
        st, ed = truvari.entry_boundaries(entry)
        if start <= st and ed <= end:
            sz = truvari.entry_size(entry)
            if sz >= 50:
                counts[1] += 1
            elif sz >= 5:
                counts[0] += 1
import json
print(json.dumps(count, indent=4))
