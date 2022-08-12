in_files=$@
mkdir -p jobs
mkdir -p utmos_files
for i in $in_files
do
    bname=$(basename $i)
    job_name=jobs/${bname%.vcf.gz}_cvt.sh
    echo "#!/bin/bash" > ${job_name}
    echo "bcftools view -i \"SVLEN == '.'\" $i | bcftools filter -S . -e \"FT == 'FAIL'\" | utmos convert --lowmem /dev/stdin utmos_files/${bname}_utmos_small.jl" >> ${job_name}
    echo "bcftools view -i \"SVLEN != '.'\" $i | bcftools filter -S . -e \"FT == 'FAIL'\" | utmos convert --lowmem /dev/stdin utmos_files/${bname}_utmos_big.jl" >> ${job_name}
done
