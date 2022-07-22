gap_bed=$1
genome=$2
bed=$3

if [ -z ${gap_bed} ]
then
    echo "no gap_bed file"
    exit
fi

if [ -z ${genome} ]
then
    echo "no genome file"
    exit
fi

if [ -z ${bed} ]
then
    echo "no bed file"
    exit
fi


bedtools slop -b 5000 -i <(zgrep -v "_" $gap_bed) -g $genome > gap_slop.bed
bedtools subtract -a $bed -b gap_slop.bed -A \
    | bedtools sort \
    | ../scripts/bed_stats.py final \
    | bgzip > tr_regions.bed.gz
tabix tr_regions.bed.gz

