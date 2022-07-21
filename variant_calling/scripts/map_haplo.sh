#!/bin/bash 

# Input fasta haplotype
fasta=$1
# Reference
ref=$2
sample_name=$3
# Directory to write the output folder
out_dir=$4

# Don't use static paths, but whatever for now
mm=/users/u233287/scratch/misc_software/minimap2-2.24/minimap2
pf=/users/u233287/scratch/misc_software/minimap2-2.24/misc/paftools.js
# These would need to be a parameter
params='-cx asm20 -m 10000 -z 10000,50 -r 50000,2000000 --end-bonus=100 --rmq=yes -O 5,56 -E 4,1 -B'
threads=8
# minimum quality score of alignments to consider
min_qual=60
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p ${out_dir}
cd ${out_dir}
anno=$DIR/annotate_cov.py
# -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 .... from pangenie
# from GIAB -z 200000,10000
# original $mm -cx asm5  -t8 -k20 --secondary=no --cs ${ref} ${fasta} | sort -k6,6 -k8,8n > ${aln_name}.paf
$mm ${params} -t${threads} --secondary=no --cs ${ref} ${fasta} \
    | sort -k6,6 -k8,8n > aln.paf
$pf stat aln.paf > aln.paf.stats.txt
cat aln.paf | $pf call -f ${ref} -q ${min_qual} -L10000 - -s ${sample_name} | vcf-sort | bgzip > aln.vcf.gz
tabix aln.vcf.gz

awk -v mq=${min_qual} '{if ($12 >= mq) print $6 "\t" $8 "\t" $9}' aln.paf > aln.bed
bedtools genomecov -i aln.bed -g ${ref}.fai -bga > cov.bed
awk '{if ($4 == 1) print $0}' cov.bed | bedtools sort > single_cov.bed

python3 $anno aln.vcf.gz cov.bed | bgzip > aln.covanno.vcf.gz
tabix aln.covanno.vcf.gz
