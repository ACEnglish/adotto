GIAB TR Benchmark Strawman
==========================

We collected TR regions for chromosome 20 from adotto_TRregions v1. (see README.md for details)
We subset these regions to those confidently covered by two alignments of the HPRC assembly.
See [adotto github](https://github.com/ACEnglish/adotto/blob/main/strawman_chr20/scripts/make_hq_bed.sh) for commands used
to create the confidently covered bed files.

The HG002 variants within regions and over 5bp were then run through `truvari anno trf`. Heuristics were applied to 
set the FILTER field to PASS or NONTR based on the annotations applied by `truvari anno trf`, variants having a 
TRFperiod > 1, and variant sequence entropy being over 0.25.
See [this script](https://github.com/ACEnglish/adotto/blob/main/strawman_chr20/scripts/tr_identification_heuristics.py) for
details.

The confidently covered regions were then split into three sets based on their intersection to HG002 variants:
* Green: A single >=5bp INDEL is present in the pVCF over the region
* Blue: More than one >=5bp INDEL is present in the pVCF over the region
* Controls: No >=5bp INDEL is present in the pVCF over the region

Each set has two associated files. 
* \[green|blue|controls\].bed.gz - The TR regions for the set
* chr20.HG002.\[green|blue|controls\].vcf.gz - The HG002 variants contained within the regions.

These sets were then consolidated to create a single benchmarking resource via
```
zcat *bed.gz | bedtools sort | bgzip > strawman.bed.gz
vcf-concat *.vcf.gz | bgzip > chr20.HG002.strawman.vcf.gz
```

To benchmark against the regions, use the development branch of version of truvari (v4.0-dev)

```
git clone https://github.com/ACEnglish/truvari.git
cd truvari
python3 -m pip install .
```

The two commands you'll need are:
```
truvari bench -b chr20.HG002.strawman.vcf.gz -c <your_calls.vcf.gz> --includebed strawman.bed.gz \
	--sizemin 5 --pick ac -o bench_result/
truvari refine --use-original --reference grch38.fasta bench_result/
```

See [truvari wiki](https://github.com/ACEnglish/truvari/wiki) for details on the outputs created.

