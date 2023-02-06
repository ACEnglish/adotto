in_vcf=$1
in_bed=$2
out_vcf=$3
REF=~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
truvari anno trf -t 4 -i $in_vcf -m 5 -f $REF \
    -T='3 7 7 80 5 5 500 -h -ngs' -r $in_bed | vcf-sort | truvari anno svinfo -m 0 | bgzip > tmp_$out_vcf
tabix tmp_$out_vcf
python annotate_tr_filter.py tmp_$out_vcf $in_bed | vcf-sort | bgzip > $out_vcf
tabix $out_vcf
rm tmp_$out_vcf
rm tmp_${out_vcf}.tbi
