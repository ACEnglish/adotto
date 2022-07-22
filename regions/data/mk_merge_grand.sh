genome=$1

if [ -z ${genome} ]
then
    echo """
Expected a ref.genome file
        The genome file should tab delimited and structured as follows:

        <chromName><TAB><chromSize>

        For example, Human (hg19):
        chr1    249250621
        chr2    243199373
        ...
        chr18_gl000207_random   4262
"""
exit
fi
echo "Initial merge"
zcat */merged.bed.gz \
    | ../scripts/bed_stats.py \
    | bedtools sort -i \
    | bedtools merge \
    | ../scripts/merged_bed_filter.py \
    | bgzip > merged.bed.gz
tabix merged.bed.gz

echo "Slop merge"
bedtools slop -i merged.bed.gz -b 25 -g $genome \
    | ../scripts/bed_stats.py merged \
    | bedtools sort \
    | bedtools merge \
    | ../scripts/merged_bed_filter.py slop \
    | bgzip > merged.slop25.bed.gz
tabix merged.slop25.bed.gz

# Get final stats real quick
zcat merged.slop25.bed.gz | ../scripts/bed_stats.py slop > /dev/null
