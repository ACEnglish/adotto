#!/bin/bash

for i in HG002 NA24385 li:NA24385
do
    echo "bcftools view -s $i GRCh38.variants.sqoff.vcf.gz | vcf-subset -e | bgzip > ${i}_comp.vcf.gz"
done | xargs -I {} -P 4 bash -c {}

for i in *comp.vcf.gz;
do
    tabix $i
done
