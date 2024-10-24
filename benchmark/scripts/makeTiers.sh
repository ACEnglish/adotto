# Get those header lines
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$DIR/../inputs/
VARDIR=$DIR/../variants/
TMPDIR=$DIR/../temp/
REPDIR=$DIR/../replicate_bench/

REF=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
PVCF=$INDIR/adotto_variants.grch38.sqoff.vcf.gz

# replaced by a manual pandas join so we can deal with new 'nulls' for some states
#python $DIR/consolidate_refine_regions.py $REPDIR/bench_ei/refine.regions.txt \
                                          #$REPDIR/bench_li/refine.regions.txt \
                                          #$REPDIR/bench_al/refine.regions.txt \
                                          #> $TMPDIR/all_refine.bed
# If you need to remake the curated states
#cut -f4,5,6 all_refine.bed | sort | uniq -c | sort -nr > states.txt

python $DIR/make_tiers.py $TMPDIR/all_refine.bed \
                          $INDIR/replicate_curated_states.txt \
                          $REF \
                          $PVCF \
                          GIABTR.HG002.benchmark.regions.bed
