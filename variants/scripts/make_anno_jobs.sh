mkdir -p jobs
mkdir -p logs
mkdir -p anno_parts
ref=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
annos=/users/u233287/scratch/code/adotto/regions/adotto_TRannotations_v0.2.bed.gz
for in_name in parts/*.bcf.gz
do  
    m_name=$(basename ${in_name})
    m_name=${m_name%.bcf.gz}
    out_name=anno_parts/${m_name}.vcf.gz
    job_name=jobs/trf_${m_name}.sh
    echo "#!/bin/bash" > $job_name
    echo "truvari anno trf -m 5 -i ${in_name} -f ${ref} -r ${annos} -t 24 | vcf-sort | bgzip > ${out_name}" >> $job_name
done
