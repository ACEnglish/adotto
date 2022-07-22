zcat grch38.simpleRepeat.truvari.bed.gz \
    | cut -f1-3 \
    | bedtools merge \
    | bedtools sort \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz
