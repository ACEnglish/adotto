
t=2
reR=/Users/english/code/regione_rust/target/release/regione_rust
for random in novl circle
do
    time $reR -g data/grch38.genome.txt \
        -A data/adotto_TRcatalog_v1.1.bed \
        -B data/grch38.protein_coding_transcript.bed \
        --mask data/grch38.exclude_regions.bed \
        -t $t -n 10000 --per-chrom --random ${random} \
        -o results/transcript_perchrom${random}1000.json 

    time $reR -g data/grch38.genome.txt \
        -A data/adotto_TRcatalog_v1.1.bed \
        -B data/grch38.epd_promoters.bed \
        --mask data/grch38.exclude_regions.bed \
        -t $t -n 10000 --per-chrom --random ${random} \
        -o results/promoters_perchrom${random}1000.json 

    time $reR -g data/grch38.genome.txt \
        -A data/adotto_TRcatalog_v1.1.bed \
        -B data/grch38.protein_coding_introns.bed \
        --mask data/grch38.exclude_regions.bed \
        -t $t -n 10000 --per-chrom --random ${random} \
        -o results/intron_perchrom${random}1000.json 

    time $reR -g data/grch38.genome.txt \
        -A data/adotto_TRcatalog_v1.1.bed \
        -B data/grch38.protein_coding_nonintrons.bed \
        --mask data/grch38.exclude_regions.bed \
        -t $t -n 10000 --per-chrom --random ${random} \
        -o results/nonintron_perchrom${random}1000.json 
done


