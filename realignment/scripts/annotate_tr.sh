#!/bin/bash

in_vcf=$1
out_name=$2
truvari anno trf -i $in_vcf \
    -e  /users/u233287/scratch/misc_software/trf409.linux64 \
    -s /users/u233287/scratch/giabvntr/bed_files/baylor/grch38.simpleRepeat.truvari.bed.gz \
    -f /users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa \
    -m 5 -t 16 | vcf-sort | bgzip > $out_name
tabix $out_name

#truvari anno hompct -i data/GRCh38.hconly.trf.vcf.gz \
#-m 5 -c 5 | vcf-sort | bgzip > data/GRCh38.hconly.trf.hom.vcf.gz
#tabix data/GRCh38.hconly.trf.hom.vcf.gz

