Strawman Process:

1 - Select the TRregions on chr20 that are covered by the HG002 assembly
	https://github.com/ACEnglish/adotto/commit/92f85c8a68bcef65a41b68224d70de8622fffb5b
	Specifically, `scripts/make_hq_bed.sh`

2 - Extract only HG002 from pVCF. (chr20)
	bcftools view -r chr20 -s HG002 -c 1

3 - Extract only chr20 from TRregions_v1.1

4 - Run truvari anno trf on the HG002 VCF
	truvari anno trf -i strawman_AllHG002.vcf.gz -o strawman_TRregions_raw.vcf -r
	adotto_TRregions_v1.1_HPRC_HG002_Covered_annotations.bed.gz -R -f
	~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa -m 5 -t 8 -T="3 7 7 80 5 5 500 -h
	-ngs"

5 - Run variants/scripts/tr_identification_heuristics.py

6 - Don't worry about green vs blue? How did I do that stuff.. I bet that's where the problem lies..?
  - I should have done that after doing the tr_identification, anyway
