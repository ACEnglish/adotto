ref=$1
vcf=$2
region=$3
# Number of bases up/down stream to add to consensus sequence (without variants)
buffer=100
#ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
#vcf=/users/u233287/scratch/insertion_ref/msru/data/inter_merge/grch38/exact/exact.vcf.gz
#vcf=/users/u233287/scratch/code/adotto/pVCFs/GRCh38.variants.squareoff.vcf.gz
#region=$1 
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "MSA realignment of region" $region
mkdir -p msa_$region

read chr st ed <<< $(echo $region | sed "s/[:-]/\t/g")
st=$(expr $st - $buffer)
ed=$(expr $ed + $buffer)
expand_region=${chr}:${st}-${ed}

echo $expand_region
# -s here
bcftools view -c 1 -r $region $vcf \
    | bcftools +fill-from-fasta /dev/stdin -- -c REF -f $ref \
    | bgzip > msa_${region}/variants.vcf.gz

tabix msa_${region}/variants.vcf.gz

python $DIR/get_reference.py $ref $expand_region > msa_${region}/haps.fa
for i in $(bcftools view -h msa_${region}/variants.vcf.gz | grep -m1 '#CHROM' | cut -f10-)
do
    samtools faidx $ref $expand_region \
        | bcftools consensus -H1 --sample $i msa_${region}/variants.vcf.gz \
        | python $DIR/fa_rename.py ${i}_1 >> msa_${region}/haps.fa
    samtools faidx $ref $expand_region \
        | bcftools consensus -H2 --sample $i msa_${region}/variants.vcf.gz \
        | python $DIR/fa_rename.py ${i}_2 >> msa_${region}/haps.fa
done
#python $DIR/remove_redundant.py msa_${region}/haps.fa > msa_${region}/haps_noredund.txt


#/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --auto msa_${region}/haps_noredund.txt > msa_${region}/aln_results.txt
/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --retree 2 --maxiterate 0 msa_${region}/haps.fa > msa_${region}/aln_results.txt
#/users/u233287/scratch/misc_software/mafft-linux64/mafft.bat --globalpair --maxiterate 1000 msa_${region}/haps_noredund.txt > msa_${region}/aln_results.txt

#./ProGraphMSA+TR.sh -o result_${region}.txt -R haps_noredund_${region}.txt

python $DIR/msa2vcf.py msa_${region}/aln_results.txt \
    | vcf-sort \
    | bcftools +fill-tags | bgzip > msa_${region}/result.vcf.gz
    #| bcftools norm -m- -c s -f $ref \
tabix msa_${region}/result.vcf.gz

echo Reports
echo "Original" >> msa_${region}/report.txt
bash $DIR/mend_report.sh msa_${region}/variants.vcf.gz >> msa_${region}/report.txt
echo "Realigned" >> msa_${region}/report.txt
bash $DIR/mend_report.sh msa_${region}/result.vcf.gz >> msa_${region}/report.txt

# Turning off validation checking for now
#python $DIR/get_reference.py $ref $region > msa_${region}/haps_final.fa
#for i in $(zgrep -m1 '#CHROM' $vcf | cut -f10-)
#do
    #samtools faidx $ref $region \
        #| bcftools consensus -H1 --sample $i msa_${region}/result.vcf.gz \
        #| python $DIR/fa_rename.py ${i}_1 >> msa_${region}/haps_final.fa
    #samtools faidx $ref $region \
        #| bcftools consensus -H2 --sample $i msa_${region}/result.vcf.gz \
        #| python $DIR/fa_rename.py ${i}_2 >> msa_${region}/haps_final.fa
#done
#
#echo "md5sums" $(sort msa_${region}/haps.fa | md5sum) $(sort msa_${region}/haps_final.fa | md5sum)


