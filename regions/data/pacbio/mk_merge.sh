bedtools sort -i repeat_catalog.hg38.bed.gz \
    | ../../scripts/bed_stats.py \
    | bedtools merge \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
