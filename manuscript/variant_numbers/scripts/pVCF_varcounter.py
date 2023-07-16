import pysam
import truvari

catalog = "/users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.1.bed.gz"
variants = "/users/u233287/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.vcf.gz"
vcf = pysam.VariantFile(variants)

fh = truvari.opt_gz_open(catalog)
print('chrom\tstart\tend\tsnp\tfive\tfifty\tsv\ttotal_len\ttotal_count')
for line in fh:
    data = line.strip().split('\t')
    chrom, start, end = data[:3]
    start = int(start)
    end = int(end)
    # snp, [1,5), [5, 50), 50+
    total_count = 0
    snp = 0
    five = 0
    fifty = 0
    sv = 0
    total_len = 0
    for entry in vcf.fetch(chrom, start, end):
        ent_start, ent_end = truvari.entry_boundaries(entry)
        if not (start <= ent_start and ent_end <= end):
            continue
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
    print(chrom, start, end, snp, five, fifty, sv, total_len, total_count, sep='\t')

