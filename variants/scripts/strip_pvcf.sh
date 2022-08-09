# Make the vcf smaller by removing/shortening some annotations and saving a compressed BCF
bcftools annotate -x FORMAT/BPDP,INFO/AF,INFO/MAF,INFO/HWE,INFO/ExcHet adotto_variants.grch38.sqoff.vcf.gz \
    | sed 's/AIL//g' | sed 's/ASS:/:/g' | bcftools view --no-version -O b -o adotto_variants.grch38.sqoff.bcf.gz
