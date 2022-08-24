#
# Given a pre-computed MSA, a comparison VCF and 
# a region, build the consensus sequence for the haplotypes of the comparison vcf
# and use mafft --add to compare the haplotypes
# Then, use a custom script to.. count matching variants by GT <-- this is where it gets weird, also
# Like, how would you count TPs/FPs etc. You could do a report for
# It was this number of variants, but after MSA it's that number of variants
# then you can just measure the similarities
# It would be nice to be able to trace the TPs/ FPs back to the representation, though.
# 

# 
# TODO - turn this into a tool
# can probably pull the region from aln_results.txt..?
# Need to be able to do parameter parsing here and in msa_realignment
# The bad part about all this is
# 1 - can't do whole genome MSA. so there's no reference
# 2 - can't keep input/output variant counts the same
# 3 - Maybe there's some way to tell that, Hey, you got FNs but we think that a more refined
# comparison will resolve them. Go use the msa...
# 
set -e
m_msa=aln_results.txt
ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
region="chr1:26399665-26399952"
buffer=100
compare_vcf=/users/u233287/scratch/code/adotto/variants/benchmarking/truthsets/HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz
compare_vcf=~/scratch/giabvntr/calls/ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/Leicester_HipSTR_GangSTR_05182022/GangSTR/Ash_trio_gangstr_100x_.vcf.gz

DIR=~/scratch/code/adotto/realignment/scripts
bcftools view -s HG002 -r $region $compare_vcf -o comp.vcf.gz -O z
tabix comp.vcf.gz

read chr st ed <<< $(echo $region | sed "s/[:-]/\t/g")
st=$(expr $st - $buffer)
ed=$(expr $ed + $buffer)
expand_region=${chr}:${st}-${ed}
echo $expand_region
samtools faidx $ref $expand_region \
    | bcftools consensus -H1 comp.vcf.gz \
    | python $DIR/fa_rename.py hap_1 > comp.fa
samtools faidx $ref $expand_region \
    | bcftools consensus -H2 comp.vcf.gz \
    | python $DIR/fa_rename.py hap_2 >> comp.fa

~/scratch/misc_software/mafft-linux64/mafft.bat --add comp.fa aln_results.txt > to_benchmark.txt

python $DIR/msa2vcf.py to_benchmark.txt \
    | vcf-sort | bgzip > bench_result.vcf.gz

bcftools query -f "[%GT ]\n" bench_result.vcf.gz | sort | uniq -c 
