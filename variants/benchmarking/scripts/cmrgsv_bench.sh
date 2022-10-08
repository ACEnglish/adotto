ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
rtg=/hgsc_software/rtg-tools/rtg-tools-3.12.1/rtg
rtg_ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.sdf

declare -A comp_vcfs
comp_vcfs[hprc]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/HG002_comp.vcf.gz
comp_vcfs[li]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/li:NA24385_comp.vcf.gz
comp_vcfs[eichler]=/users/u233287/scratch/code/adotto/variants/benchmarking/inputs/NA24385_comp.vcf.gz

bdir=/users/u233287/scratch/code/adotto/variants/benchmarking/
cmrg_vcf=$bdir/truthsets/HG002_GRCh38_CMRG_SV_single.vcf.gz
cmrg_bed=$bdir/truthsets/GRCh38_CMRG_benchmark_gene_coordinates.bed


odir=cmrg_comb_sv/
mkdir -p $odir
for samp_name in "${!comp_vcfs[@]}"
do
    m_vcf="${comp_vcfs[$samp_name]}"
    $rtg vcfeval -t $rtg_ref --squash-ploidy -e $cmrg_bed \
                 -b $cmrg_vcf -c $m_vcf -o $odir/rtg_${samp_name}

    truvari bench -b $cmrg_vcf -c $m_vcf \
                  --includebed $cmrg_bed --passonly \
                  -o $odir/truvari_d_${samp_name}
    truvari bench -b $cmrg_vcf -c $m_vcf --multimatch \
                  --includebed $cmrg_bed --passonly \
                  -o $odir/truvari_m_${samp_name}

done
