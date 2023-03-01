For a user-friendly description of how the strawman is created, see DataDescription.md

Strawman Creation process

1. Select the TRregions on chr20 that are covered by the HG002 assembly
```
	scripts/make_hq_bed.sh
```
2. Extract only HG002 from pVCF. (chr20)
```
	bcftools view -r chr20 -s HG002 -c 1
```
3. Extract only chr20 from TRregions_v1.1 (just grep)
4. Run truvari anno trf on the HG002 VCF
```
	truvari anno trf -i strawman_AllHG002.vcf.gz -o strawman_TRregions_raw.vcf \
		-r adotto_TRregions_v1.1_HPRC_HG002_Covered_annotations.bed.gz -R \
		-f ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa \
		-m 5 -t 8 -T="3 7 7 80 5 5 500 -h -ngs"
```

5. Run `scripts/tr_identification_heuristics.py` on the previously created VCF to set FILTER field
6. (Optionally) Split the regions into control/green/blue tiers using `scripts/tier_maker.py`

