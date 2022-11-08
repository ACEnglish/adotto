mkdir -p initial_alignments
# Hard coded :(
AGCFILE=/users/u233287/scratch/lrs_variants/hprc/HPRC-yr1.agc
REF=/users/u233287/scratch/insertion_ref/msru/data/reference/chm13_v2.0/chm13v2.0.fa

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p jobs
NAMES=$(~/scratch/misc_software/agc-2.1_x64-linux/agc listset ~/scratch/lrs_variants/hprc/HPRC-yr1.agc | grep -vw "CHM13Y\|HG002")
for i in $NAMES
do 
echo "#!/bin/bash" > jobs/aln_hprc_${i}.sh
bash $DIR/agc_map_haplo.sh \
    $AGCFILE \
    $i \
    $REF \
    initial_alignments/ >> jobs/aln_hprc_${i}.sh
done
