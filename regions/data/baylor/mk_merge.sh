zcat grch38.simpleRepeat.truvari.bed.gz \
    | cut -f1-3 \
    | ../../scripts/bed_stats.py \
    | bedtools merge \
    | bedtools sort \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
