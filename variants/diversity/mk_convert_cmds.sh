mkdir -p jobs
mkdir -p utmos_files
for i in ../data/parts/*.bcf.gz
do
    bname=$(basename $i)
    job_name=jobs/${bname%.bcf.gz}_cvt.sh
    echo "#!/bin/bash" > ${job_name}
    echo "bcftools view $i | utmos convert --lowmem /dev/stdin  utmos_files/${bname}_utmos_all.jl" >> ${job_name}
    echo "bcftools view -i \"SVLEN == '.'\" $i | utmos convert --lowmem /dev/stdin utmos_files/${bname}_utmos_small.jl" >> ${job_name}
    echo "bcftools view -i \"SVLEN != '.'\" $i | utmos convert --lowmem /dev/stdin utmos_files/${bname}_utmos_big.jl" >> ${job_name}
done
