hap1=$1
bed1=$2
hap2=$3
bed2=$4
sample=$5
output=$6

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
bcftools merge --force-samples -m none ${hap1} ${hap2} \
        | python ${DIR}/single_sample_cov.py /dev/stdin $bed1 $bed2 \
        | cut -f1-10 \
        | truvari anno svinfo -m 10 \
        | sed "s/sample/${sample}/" \
        | bcftools annotate -x INFO/QNAME,INFO/QSTART,INFO/QSTRAND \
        | bgzip > ${output}
tabix ${output}
