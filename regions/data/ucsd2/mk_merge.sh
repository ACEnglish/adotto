bedtools sort -i GIAB_adVNTR_short_VNTR_regions.bed.gz \
    | ../../scripts/bed_stats.py \
    | bedtools merge \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
