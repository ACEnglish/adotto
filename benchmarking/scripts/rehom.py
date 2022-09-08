"""
TRGT will report redundant alleles for homozygous variants
For example a REF ALT GT of :
A  ATAT,ATAT 1|1

This is a problem because splitting the multi allelics creates a REF/ALT of:
A ATAT 0,1
A ATAT 1,0

This script will combine these variants back into a homozygous call
"""
import sys
import pysam
import truvari
v = pysam.VariantFile(sys.argv[1])
o = pysam.VariantFile('/dev/stdout', 'w', header=v.header)
prev = next(v)
for entry in v:
    # is het redund
    if prev.pos == entry.pos and prev.chrom == entry.chrom and prev.ref == entry.ref and prev.alts[0] == entry.alts[0] \
       and truvari.get_gt(entry.samples[0]["GT"]) == truvari.GT.HET \
       and truvari.get_gt(prev.samples[0]["GT"]) == truvari.GT.HET:
        prev.samples[0]["GT"] = (1, 1)
        o.write(prev)
        # this will be an edge case problem if the last pair of variants are these split homs
        prev = next(v)
    else:
        o.write(prev)
        prev = entry
# get the last one out
o.write(prev)

