set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$DIR/../inputs/
OUTDIR=$DIR/../variants/
TMPDIR=$DIR/../temp/

mkdir -p $OUTDIR

VCF=$INDIR/adotto_variants.grch38.sqoff.vcf.gz
bcftools view -c 1 $VCF -O z -s HG002 -o $OUTDIR/HG002.vcf.gz & 
bcftools view -c 1 $VCF -O z -s NA24385 -o $OUTDIR/ei.NA24385.vcf.gz &
bcftools view -c 1 $VCF -O z -s li:NA24385 -o $OUTDIR/li.NA24385.vcf.gz &
wait

tabix $OUTDIR/HG002.vcf.gz &
tabix $OUTDIR/ei.NA24385.vcf.gz &
tabix $OUTDIR/li.NA24385.vcf.gz &
wait
