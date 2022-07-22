
zcat ./ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz \
    ./ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_lt51bp_slop5.bed.gz \
    ./ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_gt100bp_slop5.bed.gz \
    ./ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_gt10000bp_slop5.bed.gz \
    ./ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz \
    | grep -v "#" \
    | ../../scripts/bed_stats.py \
    | bedtools sort \
    | bedtools merge \
    | ../../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz

