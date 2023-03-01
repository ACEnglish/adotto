zcat /stornext/snfs4/next-gen/scratch/english/round2/code/adotto/intersection/StrawMan1.regions.bed.gz | awk '{print $1
"\t" $2 "\t" $3 "\tbench"}' > bench.bed


cat bench.bed gangstr.bed hipstr.bed trgt.bed | bedtools sort | bedtools merge -c 4 -o distinct > all_merged.bed


python make_upset.py
