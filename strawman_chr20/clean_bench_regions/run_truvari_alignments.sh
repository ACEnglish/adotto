COMP=/users/u233287/scratch/code/adotto/benchmarking/truthsets/trash/thfa_hg002.vcf.gz
BASE=chr20.HG002.replicates.vcf.gz
BED=../chr20_strawman_files_Feb.28.2023/strawman_regions.bed.gz
REF=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa

truvari bench --no-ref a --includebed $BED  -b $BASE -c $COMP --bSample HG002 -s 5 --pick ac -o bench_thfa
truvari refine --use-original --reference $REF -t 8 bench_thfa
