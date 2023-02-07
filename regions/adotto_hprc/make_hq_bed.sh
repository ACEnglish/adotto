# Select for only single coverage
genome=~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.genome_bedtools.bed
ad_mat=/users/u233287/scratch/lrs_variants/initial_alignment/giab_HG002_mat/cov.bed
ad_pat=/users/u233287/scratch/lrs_variants/initial_alignment/giab_HG002_pat/cov.bed
dipcall_hq=/users/u233287/scratch/code/adotto/regions/nist_subsetting/GRCh38_HG2-HPRC-20211005_dipcall-z2k.benchmark.bed

trr=/users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.1.bed.gz

mkdir -p temp

awk '$4 == 1' $ad_mat | cut -f1-3 > temp/hq_mat_cov.bed
awk '$4 == 1' $ad_pat | cut -f1-3 > temp/hq_pat_cov.bed

cat temp/hq_mat_cov.bed temp/hq_pat_cov.bed | bedtools sort | bedtools genomecov -i - -g $genome -bga \
    | awk '$4 == 2' | cut -f1-3 > temp/adotto_hq_cov.bed

grep -v "_" $dipcall_hq > temp/dipcall_hq_cov.bed

cat temp/adotto_hq_cov.bed temp/dipcall_hq_cov.bed | bedtools sort \
    | bedtools genomecov -i - -g $genome -bga \
    | awk '$4 == 2' | cut -f1-3 > HPRC_HG002_dipcalladotto_HQCov_spans.bed


bedtools intersect -u -f 1 -a $trr -b HPRC_HG002_dipcalladotto_HQCov_spans.bed | cut -f1-3 > temp/both_subset.bed
bedtools intersect -u -f 1 -a $trr -b temp/adotto_hq_cov.bed | cut -f1-3 > temp/adotto_subset.bed
bedtools intersect -u -f 1 -a $trr -b temp/dipcall_hq_cov.bed | cut -f1-3 > temp/dipcall_subset.bed
