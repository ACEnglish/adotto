candidates=data/tr_regions.bed.gz
mkdir -p temp

declare -A RENAMER
RENAMER['abundant_pCGG']='pCGG'
RENAMER['baylor']='UCSC'
RENAMER['giab']="GIAB"
RENAMER['pacbio']="Illumina"
RENAMER['pbsv']="pbsv"
RENAMER['trgt']="TRGT"
RENAMER['ucsd1']="UCSD1"
RENAMER['ucsd2']="UCSD2"
RENAMER['usc']="USC"

for i in data/*/merged.bed.gz
do
    src_name=$(basename $(dirname $i))
    src_name=${RENAMER[$src_name]}
    echo $src_name
    bedtools intersect -u -wa -a $candidates -b $i > temp/${src_name}.txt
done

python make_upset_data.py
