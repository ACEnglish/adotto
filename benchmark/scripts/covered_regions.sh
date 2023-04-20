# Select for confidently covered
# paths are relative to this script
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$DIR/../inputs/
OUTDIR=$DIR/../covered/
TMPDIR=$DIR/../temp/

genome=$INDIR/GRCh38_1kg_mainchrs.genome_bedtools.bed
ad_mat=$INDIR/HPRC.HG002.mat.cov.bed
ad_pat=$INDIR/HPRC.HG002.pat.cov.bed
dp_cov=$INDIR/GRCh38_HG2-HPRC-20211005_dipcall-z2k.benchmark.bed
par_reg=$INDIR/GRCh38.chrX.PAR.bed

trr=$INDIR/adotto_TRregions_v1.1.bed.gz

mkdir -p $TMPDIR
mkdir -p $OUTDIR

awk '$4 == 1' $ad_mat | cut -f1-3 > $TMPDIR/hq_mat_cov.bed
awk '$4 == 1' $ad_pat | cut -f1-3 > $TMPDIR/hq_pat_cov.bed

cat $TMPDIR/hq_mat_cov.bed $TMPDIR//hq_pat_cov.bed \
    | bedtools sort \
    | bedtools genomecov -i - -g $genome -bga > $TMPDIR/hq_diploid_cov.bed


# Non-chrX/chrY requires 1x per haplotype
awk '$4 == 2 && $1 != "chrX" && $1 != "chrY"' $TMPDIR/hq_diploid_cov.bed | cut -f1-3 > $TMPDIR/adotto_hq_cov.bed
# chrY requires 1x
awk '$4 == 1 && $1 == "chrY"' $TMPDIR/hq_diploid_cov.bed | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed

# chrX requires 2x on PAR regions
# 1) take out just chrX
grep chrX $TMPDIR/hq_diploid_cov.bed > $TMPDIR/hq_diploid_cov.chrX.bed
# 2) subtract PAR from chrX coverage and get only 1x covered
bedtools subtract -a $TMPDIR/hq_diploid_cov.chrX.bed -b $par_reg \
    | awk '$4 == 1' \
    | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed

# 3) compliment PAR to make non-PAR
bedtools complement -i $par_reg -g $genome > $TMPDIR/nonpar.bed
# 4) subtract non-PAR from chrX coverage and get only 2x covered
bedtools subtract -a $TMPDIR/hq_diploid_cov.chrX.bed -b $TMPDIR/nonpar.bed \
    | awk '$4 == 2' \
    | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed
# note, these won't be sorted, but we sort in a couple steps

# Clean out alt contig bed entries from dipcall bed
grep -v "_" $dp_cov > $TMPDIR/dipcall_hq_cov.bed

# Intersect adotto/dipcall HQ regions - if coverage == 2, they both like the region
cat $TMPDIR/adotto_hq_cov.bed $TMPDIR/dipcall_hq_cov.bed | bedtools sort \
    | bedtools genomecov -i - -g $genome -bga \
    | awk '$4 == 2' | cut -f1-3 > $OUTDIR/HPRC_HG002_covered_spans.bed

# Create TRregion subsets
bedtools intersect -u -f 1 -a $trr -b $OUTDIR/HPRC_HG002_covered_spans.bed > $OUTDIR/adotto_TRregions_v1.1_HPRC_HG002_covered.bed
# I use this with visualization
bedtools intersect -u -f 1 -a $trr -b $TMPDIR/adotto_hq_cov.bed | cut -f1-3 > $TMPDIR/adotto_trr_subset.bed
bedtools intersect -u -f 1 -a $trr -b $TMPDIR/dipcall_hq_cov.bed | cut -f1-3 > $TMPDIR/dipcall_trr_subset.bed
