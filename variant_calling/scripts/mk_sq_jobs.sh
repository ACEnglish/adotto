for i in pVCFs/parts/*.vcf.gz
do
    name=$(basename $i)
    echo "#!/bin/bash" > jobs/sqoff_${name}.sh
    echo "python scripts/annotate_pvcf_cov.py $i annotree.jl coverage.jl | bgzip > pVCFs/parts_sq/$name" >> jobs/sqoff_${name}.sh
done
