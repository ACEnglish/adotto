grep -v "_" input.bed \
    | cut -f1-3 \
    | bedtools sort \
    | ../../scripts/bed_stats.py \
    | bedtools merge \
    | bedtools sort \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
