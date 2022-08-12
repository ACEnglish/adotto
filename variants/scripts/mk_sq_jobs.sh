iDIR=pVCFs/parts/
oDIR=pVCFs/parts_sq/

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p $oDIR/
for i in $iDIR/*.vcf.gz
do
    name=$(basename $i)
    echo "#!/bin/bash" > jobs/sqoff_${name}.sh
    echo "python $DIR/annotate_pvcf_cov.py $i annotree.jl coverage.jl | bcftools +fill-tags | bgzip > $oDIR/$name" >> jobs/sqoff_${name}.sh
done
