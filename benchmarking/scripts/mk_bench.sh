set -e

odir=$1
bed=$2
mkdir $odir

if [[ ! -z "$bed" ]]
then
    bed="--includebed $bed"
fi

truvari bench -s 5 -c inputs/HG002_GangSTR_100x.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o  $odir/adotto_gangstr/ $bed
truvari bench -s 5 -c inputs/HG002_HipSTR_100x.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o $odir/adotto_hipstr/ $bed
truvari bench -s 5 -c inputs/HG002_TRGT.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o $odir/adotto_trgt/ $bed

truvari bench -s 5 -c inputs/HG002_GangSTR_100x.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o  $odir/thfa_gangstr/ $bed
truvari bench -s 5 -c inputs/HG002_HipSTR_100x.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o $odir/thfa_hipstr/ $bed
truvari bench -s 5 -c inputs/HG002_TRGT.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o $odir/thfa_trgt/ $bed



