in_bed=$1
in_vcf=$2

truvari anno svinfo -i $in_vcf -m 5 | bcftools view -i "SVLEN != '.'" | bgzip > tmp.vcf.gz
bedtools intersect -c -a $in_bed -b tmp.vcf.gz | cut -f1,2,3,5
