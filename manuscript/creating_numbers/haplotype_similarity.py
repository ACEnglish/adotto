"""
Compares HG002 haplotypes to non-HG002 haplotypes from other samples in the pVCF
for each Tier1 >=5bp variant containing region in the benchmark
then run
import pandas as pd
d = pd.read_csv('tmp.txt', sep='\t', names=['chrom', 'start', 'end', 's1', 's2'])
select = (d['s1'] < 0.999) | (d['s2'] < 0.999)
cnt = select.sum()
print(cnt, cnt / len(d))
"""
import truvari
import pysam
import sys
from truvari.phab import (
    fasta_reader,
    extract_reference,
    collect_haplotypes,
    consolidate_haplotypes_with_reference
)

do_chrom = sys.argv[1]
regions = "/users/u233287/scratch/code/adotto/benchmark/GIABTR.HG002.benchmark.regions.bed.gz"
reference = "/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa"
variants = "/users/u233287/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.vcf.gz"

vcf = pysam.VariantFile(variants)
hap_jobs = []
for sample in vcf.header.samples:
    if sample not in ['NA24385', 'li:NA24385']:
        for hap in [1, 2]:
            hap_jobs.append([variants, sample, False, hap])
 
hold_fn = f"{do_chrom}_regions.txt"
with open(hold_fn, 'w') as reg_fout:
    fh = truvari.opt_gz_open(regions)
    for line in fh:
        data = line.strip().split('\t')
        # Tier1 with HG002 >=5bp 
        if data[3] == "Tier1" and int(data[5]) & 0x1 and data[0] == do_chrom:
            reg_fout.write(f"{data[0]}:{data[1]}-{data[2]}\n")

ref_haps_fn = extract_reference(hold_fn, reference)

samp_haps = collect_haplotypes(ref_haps_fn, hap_jobs, 8)
haplotypes = consolidate_haplotypes_with_reference(samp_haps, ref_haps_fn)
for seq_bytes in haplotypes:
    fasta = {k:v.decode() for k,v in fasta_reader(seq_bytes.decode(), name_entries=False)}
    ref_key = [_ for _ in fasta.keys() if _.startswith("ref_")][0]
    h1_seq = fasta[[_ for _ in fasta.keys() if _.startswith("HG002_1")][0]]
    h2_seq = fasta[[_ for _ in fasta.keys() if _.startswith("HG002_2")][0]]
    h1_sim = 0
    h2_sim = 0
    other_haps = [_ for _ in fasta.keys() if not (_.startswith("ref_") or _.startswith("HG002"))]
    for key in other_haps:
        h1_sim = max(h1_sim, truvari.seqsim(h1_seq, fasta[key]))
        h2_sim = max(h1_sim, truvari.seqsim(h2_seq, fasta[key]))
    reg = ref_key.split('_')[1].replace(':', '\t').replace('-', '\t')
    print(reg, h1_sim, h2_sim, sep='\t')
