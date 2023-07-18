import sys
# "bedtools intersect -a <( zgrep common dbSnp153_chr1.bed.gz) -b /users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.1.bed.gz -wo"
dREF=4
dALTs=6
CLS=13
NOTES=14
chrom=17
start=18
end=19

snp = 0
five = 0
fifty = 0
sv = 0
total = 0

print('chrom\tstart\tend\tsnp\tfive\tfifty\tsv\ttotal_len')
for line in sys.stdin:
    data = line.strip().split('\t')
    alt_len = [len(_) for _ in data[dALTs].split(',')][:-1]
    if not alt_len:
        continue
    max_alt_len = max(alt_len)
    sz = abs(max_alt_len - len(line[dREF]))
    # So, use ref/alts to get the variant size
    # and do it for notes to get ... rare/common (however many types there are)
    snp = 0
    five = 0
    fifty = 0
    sv = 0
    total = 0

    if sz == 0:
        snp += 1
    elif sz < 5:
        five += 1
    elif sz < 50:
        fifty += 1
    else:
        sv += 1
    total += sz
    print(data[chrom], data[start], data[end], snp, five, fifty, sv, total, sep='\t')
    

