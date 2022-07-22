bedtools sort -i ensembleTR_loci_list.txt \
    | ../../scripts/bed_stats.py \
    | bedtools merge \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
