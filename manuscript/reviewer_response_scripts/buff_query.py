import sys
import pysam
import joblib
import truvari
import numpy as np
import pandas as pd
from collections import Counter

catalog = "../regions/adotto_TRregions_v1.1.bed"
benchmark = "../benchmark/GIABTR_benchmark_v1.0/GIABTR.HG002.benchmark.regions.bed.gz"
vcf_fn = "../benchmark/GIABTR_benchmark_v1.0/GIABTR.HG002.benchmark.trfinfo.vcf.gz"
catalog = pd.read_csv(catalog, sep='\t').set_index(['chr', 'start', 'end'])
h_benchmark = pd.read_csv(benchmark, sep='\t', names=['chr', 'start', 'end', 'tier', 'repl', 'vflag', 'entropy', 'ad1',
                                                    'ad2']).set_index(['chr', 'start', 'end'])
benchmark = catalog[catalog.index.isin(h_benchmark.index)]
benchmark['my_tier'] = h_benchmark['tier']
def within(s, e, ub, db):
    return ub <= s <= e <= db

n_buff = 0
n_buff_trf = 0
n_inside = 0
n_inside_trf = 0
n_entirely_buff = 0
n_entirely_buff_trf = 0
n_entirely_buff_reg = 0
n_entirely_buff_butalso_inside = 0
vcf = pysam.VariantFile(vcf_fn)
out_vcf = pysam.VariantFile("interesting.vcf", 'w', header=vcf.header)
tier_cnt = Counter()
for _, line in benchmark.reset_index().iterrows():
    chrom, start, end = line['chr'], line['start'], line['end']
    up_buff = line['up_buff']
    dn_buff = line['dn_buff']
    has_buff = False
    has_inside = False
    for entry in vcf.fetch(chrom, start, end):
        sz = truvari.entry_size(entry)
        if sz < 5:
            continue
        s,e = truvari.entry_boundaries(entry)
        if not (start <= s <= e <= end):
            continue
        n_inside += 1
        is_tr = "TRFrepeat" in entry.info
        if is_tr:
            n_inside_trf += 1
        if not (start + up_buff <= s <= e <= end - dn_buff):
            n_buff += 1
            if is_tr:
                n_buff_trf += 1
        else:
            has_inside = True
        if within(s, e, start, start+up_buff - 1) or within(s, e, end-dn_buff+1, end):
            n_entirely_buff += 1
            has_buff = True
            if is_tr:
                n_entirely_buff_trf += 1
                out_vcf.write(entry)
    if has_buff:
        n_entirely_buff_reg += 1
        tier_cnt[line['my_tier']] += 1
        if has_inside:
            n_entirely_buff_butalso_inside += 1

print('inside:', n_inside)
print('inside_trf:', n_inside_trf)
print('buff:', n_buff)
print('buff_trf:', n_buff_trf)
print('n_entirely_buff', n_entirely_buff)
print('n_entirely_buff_trf', n_entirely_buff_trf)
print('n_entirely_buff_reg', n_entirely_buff_reg)
print("n_entirely_buff_butalso_inside", n_entirely_buff_butalso_inside)
print('tier', tier_cnt)
