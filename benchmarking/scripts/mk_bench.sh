rm -rf bench
mkdir -p bench

truvari bench -s 5 -c inputs/HG002_GangSTR_100x.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o  bench/adotto_gangstr/
truvari bench -s 5 -c inputs/HG002_HipSTR_100x.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o bench/adotto_hipstr/
truvari bench -s 5 -c inputs/HG002_TRGT.vcf.gz --passonly -b truthsets/adotto_trfanno_hg002.vcf.gz \
    -o bench/adotto_trgt/

truvari bench -s 5 -c inputs/HG002_GangSTR_100x.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o  bench/thfa_gangstr/
truvari bench -s 5 -c inputs/HG002_HipSTR_100x.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o bench/thfa_hipstr/
truvari bench -s 5 -c inputs/HG002_TRGT.vcf.gz --passonly -b truthsets/thfa_trfanno_hg002.vcf.gz \
    -o bench/thfa_trgt/



