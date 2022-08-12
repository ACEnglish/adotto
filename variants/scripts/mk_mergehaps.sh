pairs=$1
mkdir -p jobs
mkdir -p logs

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cat $pairs | while read h1 c1 h2 c2 proj samp
do
    echo "#!/bin/bash" > jobs/hapmerge_${proj}_${samp}.sh
    echo "bash $DIR/merge_haps.sh $h1 $c1 $h2 $c2 $samp hapo_merged/${proj}_${samp}.vcf.gz" >> jobs/hapmerge_${proj}_${samp}.sh
done


