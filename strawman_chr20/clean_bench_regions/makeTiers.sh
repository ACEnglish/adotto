paste bench_eichler/refine.regions.txt bench_li/refine.regions.txt \
        | cut -f1,2,3,13,26 \
        | grep FP \
        | grep -v TN > eichli_FP_regions.bed
grep "FP\|FN" bench_thfa/refine.regions.txt | cut -f1-3,13 > thfa_FNP_regions.bed
cat thfa_FNP_regions.bed eichli_FP_regions.bed > chr20_strawman_Tier2_regions.bed
bedtools intersect -v -a ../chr20_strawman_files_Feb.28.2023/strawman_regions.bed.gz \
        -b <(cut -f1-3 chr20_strawman_Tier2_regions.bed) > chr20_strawman_Tier1_regions.bed
