input=$1
output=$2
set -e
ref=~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
echo 'a'
bcftools +fill-from-fasta $input -- -c REF -f ${ref} > a
echo 'b'
bcftools norm -cs -m-any -f ${ref} a > b
echo 'e'
bcftools view b -c 1 -O z -o $output 
