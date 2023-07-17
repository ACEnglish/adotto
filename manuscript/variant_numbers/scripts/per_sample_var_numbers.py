import sys
import truvari
import joblib
import pysam
import numpy as np

do_chrom = sys.argv[1]
def bed_parser(line):
    if line == "":
        return None
    data = line.strip().split('\t')
    start = int(data[1])
    end = int(data[2])
    samps = data[4].split(',')
    is_tr = 1 if "TR" in samps else 0
    if is_tr:
        samps.remove("TR")
    return {'chrom': data[0], 'start': start, 'end': end, 'tr': is_tr, 'samps': samps}

variants = {} # key of TOTAL/TR/each-sample values of snp/five/fifty/sv/len
with open("/stornext/snfs4/next-gen/scratch/english/round2/code/adotto/manuscript/variant_numbers/all_names.txt", 'r') as fh:
    names = fh.readline().strip().split(' ')
    for sample in names[1:]:
        variants[sample] = np.zeros((2, 5), dtype=int)

# each sample: 5x2
def entry_to_cats(entry, cat_line):
    """
    """
    sz = truvari.entry_size(entry)
    idx = 0
    if sz == 0:
        idx = 0
    elif sz < 5:
        idx = 1
    elif sz < 50:
        idx = 2
    else:
        idx = 3
    is_tr = cat_line['tr']
    for sample in cat_line['samps']:
        if 1 in entry.samples[sample]['GT']:
            variants[sample][is_tr, idx] += 1
            variants[sample][is_tr, 4] += sz


vcf = pysam.VariantFile("/users/u233287/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.vcf.gz")
vcf = vcf.fetch(do_chrom)

cat = open("giant_coverage_allsamples.bed", 'r')
cur_cat = bed_parser(cat.readline())
while cur_cat['chrom'] != do_chrom:
    cur_cat = bed_parser(cat.readline())

cur_var = next(vcf)

while cur_cat is not None and cur_cat['chrom'] == do_chrom and cur_var is not None and cur_var.chrom == do_chrom:
    if cur_cat['end'] < cur_var.start:
        cur_cat = bed_parser(cat.readline())
        continue
    if cur_var.start < cur_cat['start']:
        try:
            cur_var = next(vcf)
        except StopIteration:
            cut_var = None
        continue
    # This variant is covered by an informative cat(bed) line
    entry_to_cats(cur_var, cur_cat)
    # now I've handled this variant, I need to move on to the next one
    try:
        cur_var = next(vcf)
    except StopIteration:
        cur_var = None

joblib.dump(variants, f"persamp_variants_{do_chrom}.jl")
