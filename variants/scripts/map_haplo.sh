#!/bin/bash 


# Input fasta haplotype
fasta=$(realpath $1)
# Reference
ref=$(realpath $2)
sample_name=$3
# Directory to write the output folder
out_dir=$(realpath $4)

# These would need to be a parameter / config
params='-cx asm20 -m 10000 -z 10000,50 -r 50000,2000000 --end-bonus=100 --rmq=yes -O 5,56 -E 4,1 -B'
threads=8
# minimum quality score of alignments to consider
min_qual=60

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if ! command -v minimap2 &> /dev/null
then
    echo "Program 'minimap2' not found in environment"
    exit 1
fi

if ! command -v paftools.js &> /dev/null
then
    echo "Program 'paftools.js' not found in environment"
    exit 1
fi

mkdir -p ${out_dir}
cd ${out_dir}
anno=$DIR/annotate_cov.py
minimap2 ${params} -t${threads} --secondary=no --cs ${ref} ${fasta} \
    | sort -k6,6 -k8,8n > aln.paf
paftools.js stat aln.paf > aln.paf.stats.txt
cat aln.paf | paftools.js call -f ${ref} -q ${min_qual} -L10000 -s ${sample_name} - \
    | vcf-sort \
    | bcftools +fill-from-fasta /dev/stdin -- -c REF -f ${ref} \
    | bgzip > aln.vcf.gz
tabix aln.vcf.gz

awk -v mq=${min_qual} '{if ($12 >= mq) print $6 "\t" $8 "\t" $9}' aln.paf > aln.bed
bedtools genomecov -i aln.bed -g ${ref}.fai -bga > cov.bed
awk '{if ($4 == 1) print $0}' cov.bed | bedtools sort > single_cov.bed

python3 $anno aln.vcf.gz cov.bed | bgzip > aln.covanno.vcf.gz
tabix aln.covanno.vcf.gz
