
t=4
reR=/Users/english/code/regioners/target/release/regioners
for random in novl circle
do
    #time $reR -g data/grch38.genome.txt \
        #-A data/adotto_TRcatalog_v1.1.bed \
        #-B data/grch38.protein_coding_transcript.bed \
        #--mask data/grch38.exclude_regions.bed \
        #-t $t -n 10000 --per-chrom --random ${random} \
        #-o results/transcript_perchrom${random}1000.json 

    #time $reR -g data/grch38.genome.txt \
        #-A data/adotto_TRcatalog_v1.1.bed \
        #-B data/grch38.epd_promoters.bed \
        #--mask data/grch38.exclude_regions.bed \
        #-t $t -n 10000 --per-chrom --random ${random} \
        #-o results/promoters_perchrom${random}1000.json 

    #time $reR -g data/grch38.genome.txt \
        #-A data/adotto_TRcatalog_v1.1.bed \
        #-B data/grch38.protein_coding_introns.bed \
        #--mask data/grch38.exclude_regions.bed \
        #-t $t -n 10000 --per-chrom --random ${random} \
        #-o results/intron_perchrom${random}1000.json 

    #time $reR -g data/grch38.genome.txt \
        #-A data/adotto_TRcatalog_v1.1.bed \
        #-B data/grch38.protein_coding_nonintrons.bed \
        #--mask data/grch38.exclude_regions.bed \
        #-t $t -n 10000 --per-chrom --random ${random} \
        #-o results/nonintron_perchrom${random}1000.json 

    #time $reR -g data/grch38.genome.txt \
        #-A data/adotto_TRcatalog_v1.1.bed \
        #-B data/grch38.protein_coding_CDS.bed \
        #--mask data/grch38.exclude_regions.bed \
        #-t $t -n 10000 --per-chrom --random ${random} \
        #-o results/CDS_perchrom${random}1000.json 

    time $reR -g data/grch38.genome.txt \
        -A data/adotto_TRcatalog_v1.1.bed \
        -B data/grch38.protein_coding_UTR.bed \
        --mask data/grch38.exclude_regions.bed \
        -t $t -n 10000 --per-chrom --random ${random} \
        -o results/UTR_perchrom${random}1000.json 

done


