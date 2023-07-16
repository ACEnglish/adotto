bedtools intersect -a <( zgrep common dbSnp153_chr1.bed.gz) \
    -b /users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.1.bed.gz -wo \
    | python scripts/dbsnp_varcounter.py \
    | gzip > TRcatalog_chr1_dbSNPcounts.txt.gz
