ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
rtg=/hgsc_software/rtg-tools/rtg-tools-3.12.1/rtg
rtg_ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.sdf

declare -A comp_vcfs
comp_vcfs[hprc]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/HG002_comp.vcf.gz
#comp_vcfs[li]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/li:NA24385_comp.vcf.gz
#comp_vcfs[eichler]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/NA24385_comp.vcf.gz

bdir=/users/u233287/scratch/code/adotto/variants/benchmarking/
thfa_vcf=$bdir/truthsets/HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz
thfa_bed=$bdir/truthsets/HPRC-cur.20211005-align2-GRCh38.dip.bed
thfa_sm_bed=$bdir/truthsets/temp_hc_regions/GRCh38_HG2-HPRC-20211005_dipcall-z2k.smallvar.benchmark.bed
thfa_lg_bed=$bdir/truthsets/temp_hc_regions/GRCh38_HG2-HPRC-20211005_dipcall-z2k.SV.benchmark.bed
cmrg_vcf=$bdir/truthsets/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
cmrg_bed=$bdir/truthsets/HG002_GRCh38_CMRG_SV_v1.00.bed

thfa_ma_vcf=$bdir/truthsets/HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz
cmrg_sm_vcf=$bdir/truthsets/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz
cmrg_sm_bed=$bdir/truthsets/HG002_GRCh38_CMRG_smallvar_v1.00.bed

mkdir -p no_bed/
for samp_name in "${!comp_vcfs[@]}"
do
    m_vcf="${comp_vcfs[$samp_name]}"
    #$rtg vcfeval -t $rtg_ref --squash-ploidy \
    #             -b $thfa_vcf -c $m_vcf -o no_bed/rtg_thfa_${samp_name}

    truvari bench -b $thfa_vcf -c $m_vcf --sizemin 5 --sizefilt 5 \
                  --passonly -f $ref -o no_bed/ref_single_truvari_thfa_${samp_name}
done
