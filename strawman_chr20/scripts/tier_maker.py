"""
Given the bedfile and the vcf, create the control, green, and blue tiers' bed file.
control: - no Variants >=5bp
green: - 1 variant >= 5bp
blue: - >1 variant >= 5bp
"""
import sys
import pysam
import truvari
import numpy as np

vcf = pysam.VariantFile(sys.argv[1])
tree, num_regions = truvari.build_anno_tree(sys.argv[2])

cnts = np.zeros(num_regions, dtype=int)

for entry in vcf:
    if truvari.entry_size(entry) < 5:
        continue
    hits = tree[entry.chrom].at(entry.start)
    if len(hits) == 0:
        print(entry)
        print('No hit')
        sys.exit(1)
    if len(hits) > 1:
        print(entry)
        print("Multi Hits")
        sys.exit(1)
    j = list(hits)[0]
    cnts[j.data] += 1

outs = [open("test_orig_control.bed", 'w'),
        open("test_orig_green.bed", 'w'),
        open("test_orig_blue.bed", 'w')]

for pos, region in enumerate(truvari.opt_gz_open(sys.argv[2])):
    outs[min(cnts[pos], 2)].write(region)

for i in outs:
    i.close()
