mkdir -p jobs
mkdir -p logs

cat metadata/hap_pairs.txt | while read h1 c1 h2 c2 proj samp
do
    echo "#!/bin/bash" > jobs/hapmerge_${proj}_${samp}.sh
    echo "bash scripts/merge_haps.sh $h1 $c1 $h2 $c2 $samp hapo_merged/${proj}_${samp}.vcf.gz" >> jobs/hapmerge_${proj}_${samp}.sh
done


