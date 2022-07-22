#ref=$1
#vcf=$2
#region=$3
ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
vcf=/users/u233287/scratch/insertion_ref/msru/data/inter_merge/grch38/exact/exact.vcf.gz
vcf=/users/u233287/scratch/code/adotto/pVCFs/GRCh38.variants.squareoff.vcf.gz
region=$1 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "MSA realignment of region" $region
mkdir -p results/msa_$region

#bcftools norm -m - -c s -f $ref -r $region $vcf > results/msa_${region}/variants.vcf
bcftools view -r $region -o results/msa_${region}/variants.vcf.gz -O z $vcf
tabix results/msa_${region}/variants.vcf.gz

python $DIR/get_reference.py $ref $region > results/msa_${region}/haps.fa
for i in $(zgrep -m1 '#CHROM' $vcf | cut -f10-)
do
    samtools faidx $ref $region | bcftools consensus -H1 --sample $i $vcf | python $DIR/fa_rename.py ${i}_1 >> results/msa_${region}/haps.fa
    samtools faidx $ref $region | bcftools consensus -H2 --sample $i $vcf | python $DIR/fa_rename.py ${i}_2 >> results/msa_${region}/haps.fa
done
python $DIR/remove_redundant.py results/msa_${region}/haps.fa > results/msa_${region}/haps_noredund.txt


#/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --auto results/msa_${region}/haps_noredund.txt > results/msa_${region}/aln_results.txt
/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --auto results/msa_${region}/haps.fa > results/msa_${region}/aln_results.txt
#/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --globalpair --maxiterate 1000 msa_${region}/haps_noredund.txt > msa_${region}/aln_results.txt

#./ProGraphMSA+TR.sh -o result_${region}.txt -R haps_noredund_${region}.txt

python $DIR/msa2vcf.py results/msa_${region}/aln_results.txt \
    | vcf-sort \
    | bcftools +fill-tags | bgzip > results/msa_${region}/result.vcf.gz
    #| bcftools norm -m- -c s -f $ref \
tabix results/msa_${region}/result.vcf.gz

echo Reports
bash $DIR/mend_report.sh results/msa_${region}/variants.vcf.gz
bash $DIR/mend_report.sh results/msa_${region}/result.vcf.gz

python $DIR/get_reference.py $ref $region > results/msa_${region}/haps_final.fa
for i in $(zgrep -m1 '#CHROM' $vcf | cut -f10-)
do
    samtools faidx $ref $region \
        | bcftools consensus -H1 --sample $i results/msa_${region}/result.vcf.gz \
        | python $DIR/fa_rename.py ${i}_1 >> results/msa_${region}/haps_final.fa
    samtools faidx $ref $region \
        | bcftools consensus -H2 --sample $i results/msa_${region}/result.vcf.gz \
        | python $DIR/fa_rename.py ${i}_2 >> results/msa_${region}/haps_final.fa
done

echo "md5sums" $(sort results/msa_${region}/haps.fa | md5sum) $(sort results/msa_${region}/haps_final.fa | md5sum)


