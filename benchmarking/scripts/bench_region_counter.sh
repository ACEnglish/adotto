in_dir=$1
in_bed=$2

bedtools intersect -c -a $in_bed -b $in_dir/tp-base.vcf.gz > a
bedtools intersect -c -a $in_bed -b $in_dir/tp-call.vcf.gz > b
bedtools intersect -c -a $in_bed -b $in_dir/fn.vcf.gz > c
bedtools intersect -c -a $in_bed -b $in_dir/fp.vcf.gz > d

paste a b c d | cut -f1,2,3,5,10,15,20 > $in_dir/region_counts.bed
rm a b c d
