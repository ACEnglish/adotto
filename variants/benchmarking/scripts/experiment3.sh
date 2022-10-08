ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
rtg=/hgsc_software/rtg-tools/rtg-tools-3.12.1/rtg
rtg_ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.sdf

bdir=/users/u233287/scratch/code/adotto/variants/benchmarking/
base=$bdir/cmrg_comb/phab_hprc/cmrg_giab.vcf.gz
comp=$bdir/cmrg_comb/phab_hprc/cmrg_adotto.vcf.gz
cmrg_bed=$bdir/truthsets/HG002_GRCh38_CMRG_v1.00.bed

mkdir -p cmrg_region_phab
#truvari bench -b $base -c $comp -s 0 --includebed $cmrg_bed -o cmrg_region_phab/truvari_hprc

$rtg vcfeval -t $rtg_ref --squash-ploidy -e $cmrg_bed \
             -b $base -c $comp -o cmrg_region_phab/rtg_hprc

bcftools view -R truthsets/HG002_GRCh38_CMRG_v1.00.bed cmrg_comb/phab_hprc/phab_result.vcf.gz | python scripts/phab_gtcounter.py
