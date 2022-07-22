"""
Given an input VCF
Find blocks of Variants that are isolated by at least BUFFER base pairs
"""
import sys
import pysam
import truvari

buff = 200

for chrom, start, end in truvari.vcf_ranges(sys.argv[1], buff):
    print(f"{chrom}\t{start - 200}\t{end + 200}")
