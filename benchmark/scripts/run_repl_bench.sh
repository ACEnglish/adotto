set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$DIR/../inputs/
VARDIR=$DIR/../variants/
TMPDIR=$DIR/../temp/
OUTDIR=$DIR/../replicate_bench/

mkdir -p $OUTDIR

BASEVCF=$VARDIR/HG002.vcf.gz
REP1=$VARDIR/ei.NA24385.vcf.gz
REP2=$VARDIR/li.NA24385.vcf.gz
REPA=$INDIR/HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz

REF=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
BED=$INDIR/../covered/adotto_TRregions_v1.1_HPRC_HG002_covered.bed

truvari bench --includebed $BED  -b $BASEVCF -c $REP1 -s 5 --pick ac -o $OUTDIR/bench_ei &
truvari bench --includebed $BED  -b $BASEVCF -c $REP2 -s 5 --pick ac -o $OUTDIR/bench_li &
truvari bench --includebed $BED  -b $BASEVCF -c $REPA -s 5 --pick ac -o $OUTDIR/bench_al &
wait

truvari refine --use-original --reference $REF -t 8 $OUTDIR/bench_ei
truvari refine --use-original --reference $REF -t 8 $OUTDIR/bench_li
truvari refine --use-original --reference $REF -t 8 $OUTDIR/bench_al
wait

