
VCF=chr20.HG002.replicates.vcf.gz
BED=../chr20_strawman_files_Feb.28.2023/strawman_regions.bed.gz
REF=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
truvari bench --no-ref a --includebed $BED  -b $VCF -c $VCF --bSample HG002 --cSample NA24385 -s 5 --pick ac -o bench_eichler
truvari bench --no-ref a --includebed $BED  -b $VCF -c $VCF --bSample HG002 --cSample li:NA24385 -s 5 --pick ac -o bench_li

truvari refine --use-original --reference $REF -t 8 bench_eichler
truvari refine --use-original --reference $REF -t 8 bench_li
