import pysam
import truvari

variants = "/users/u233287/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.vcf.gz"
vcf = pysam.VariantFile(variants)

print('chrom\tsnp\tfive\tfifty\tsv\ttotal_len\ttotal_count')
total_count = 0
snp = 0
five = 0
fifty = 0
sv = 0
total_len = 0

cur_chrom = None
for entry in vcf:
    if cur_chrom != entry.chrom:
        if cur_chrom is None:
            cur_chrom = entry.chrom
        else:
            print(cur_chrom, snp, five, fifty, sv, total_len, total_count, sep='\t')
            cur_chrom = entry.chrom
            total_count = 0
            snp = 0
            five = 0
            fifty = 0
            sv = 0
            total_len = 0
    total_count += 1
    sz = truvari.entry_size(entry)
    if sz == 0:
        snp += 1
    elif sz < 5:
        five += 1
    elif sz < 50:
        fifty += 1
    else:
        sv += 1
    total_len += sz
print(cur_chrom, snp, five, fifty, sv, total_len, total_count, sep='\t')

