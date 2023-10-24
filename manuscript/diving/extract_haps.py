"""
For each region, 
For each sample,
Extract the haplotype
Save as a fasta for input to trviz
"""
import sys
import joblib
import pysam
from truvari.phab import extract_reference, collect_haplotypes
ref_fn = "/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa"
bed_fn = "m_bed.bed" # These are the regions we want to visualize
vcf_fn = "/users/u233287/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.vcf.gz"
region_fn = "tmp.bed"

regions = []
with open(bed_fn) as fh, open(region_fn, 'w') as fout:
    for line in fh:
        data = line.strip().split('\t')
        fout.write("%s:%s-%s\n" % tuple(data))
        data[1] = int(data[1])
        data[2] = int(data[2])
        regions.append(data)

ref_haps_fn = extract_reference(region_fn, ref_fn)
print("# Copy this so you can reuse it")
print(ref_haps_fn)

vcf = pysam.VariantFile(vcf_fn)
jobs = []
for sample in vcf.header.samples:
    for hap in [1, 2]:
        jobs.append([vcf_fn, sample, False, hap])

sample_haplotypes = collect_haplotypes(ref_haps_fn, jobs, 4)

for k in sample_haplotypes.keys():
    sample_haplotypes[k].seek(0)
    sample_haplotypes[k] = sample_haplotypes[k].read()
joblib.dump(sample_haplotypes, 'haps.jl')

