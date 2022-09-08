In Progress

Plan
====

* Each Program’s Variant Counts / Sizes
* Advanced bench
  * can we use phab to make the performance match between truthsets
  * can we use phab to make the discovery callsets performance equal between truthsets

Why using both?
because they don’t match up directly very well meaning they use different variant representations. Adotto parameters
were tuned to maximize GIAB v0.6 Tier1 SV performance. I expect software has been tuning to this old SV standard and so
they’ve likely made their variants look more like that. Ho - no difference ,Ha - adotto matches better from the start

Why do phab?
Because the Adotto/THFA HG002 are the same sample, they only have variant representation differences. So if phab works,
then they should look the same out of the other end.

`tables/`
========
Summary data for analysis and plotting will be placed in the tables directory.
* per_region_varcounts.txt
*	Region	AdottoCount	THFACount

Collect comparison variants
===========================
GangSTR/HipSTR
--------------
```
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/Leicester_HipSTR_GangSTR_05182022/ 
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/analysis/Leicester_HipSTR_GangSTR_05182022/ 
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/analysis/Leicester_HipSTR_GangSTR_05182022/ 
```
I'm only using the AshkenazimTrio 100x to start and placing them into
```
inputs/AJTrio_GangSTR_100x.vcf.gz
inputs/AJTrio_HipSTR_100x.vcf.gz
```

These VCFs contain trios and need the variant representations changed to truvari's expectations. For example,
the line
```
#before
chr1    40620   .       ttcttcttcttc    ttcttcttc ...
#after
chr1    40619   .       TTTC    T ...
```

This is done with:
```bash
ref=~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa
bcftools +fill-from-fasta inputs/AJTrio_GangSTR_100x.vcf.gz -- -c REF -f ${ref} \
	| bcftools norm -m-any -f ${ref} \
	| bcftools view -c 1 -s HG002 -O z -o inputs/HG002_GangSTR_100x.vcf.gz
```
Note the HipSTR had to have a site at chr10:124120653 manually removed. 

TRGT
----
provided by Egor Dolzhenko @ PacBio
```
inputs/TRGT-HG002.vcf.gz
```
Same 'cleaning' as above was performed but with an additional parameter of `bcftools norm -cs` due to `Duplicate alleles at chr1:16683; run with -cw to turn the error into warning or with -cs to fix.`
Also, had to run `scripts/rehom.py` to merge redundant hets that comprise a hom. See the script for details

Setting up truthset
===================

Subset to the TRregions
-----------------------
```bash
bcftools view -R ~/scratch/code/adotto/regions/adotto_TRannotations_v0.2.bed.gz -c 1 -s HG002
	~/scratch/code/adotto/variants/data/adotto_variants.grch38.sqoff.bcf.gz \
	| bcftools filter -e "FT != 'P'" \
	| bgzip > truthsets/adotto_hg002.vcf.gz
```

And for the GIAB TrioHifiAsm (thfa)
```bash
cat GRCh38_HG2-HPRC-20211005_dipcall-z2k.*bed \
	| bedtools sort | bedtools merge > truthsets/thfa_hc.bed
bcftools view -R truthsets/thfa_hc.bed HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz \
	| bcftools view -T adotto_TRannotations_v0.2.bed.gz \
	| bgzip > truthsets/thfa_hg002.vcf.gz
```
Where here I'm subsetting to GIAB's defined 'HC' regions and then the Adotto TR regions.

Annotate truthsets
------------------
```bash
truvari anno trf -t 12 -i truthsets/adotto_hg002.vcf.gz \
	-f GRCh38_1kg_mainchrs.fa \
	-r adotto_TRannotations_v0.2.bed.gz -m 0 \
  | vcf-sort | bgzip > truthsets/adotto_trfanno_hg002.vcf.gz
```


Compare the truthsets
---------------------

TRF annotations as well as general type/size summaries.

First, let's make quick size/type summaries.

Then, we can look at the bench intersection of the truthsets

How many variants are annotated? By Type. Anything about Distribution? This is kinda difficult, actually.

And we can look at how many FN/FP there are per-region. If we see them piling up, that should be evidence that the
calls aren't matching great due to different representations.

```bash
bedtools intersect -c -a adotto_TRannotations_v0.2.bed.gz -b fn.vcf.gz > intersect_fn.txt
bedtools intersect -c -a adotto_TRannotations_v0.2.bed.gz -b fp.vcf.gz > intersect_fp.txt
```

Basic Truvari bench
-------------------
run
```
bash scripts/mk_bench.sh
```
