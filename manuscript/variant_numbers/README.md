For every Adotto_v1.1 catalog region:

Count the number of variants in the pVCF.
Record sizes and allele deltas

`scripts/pVCF_varcounter.py`

Count the number of variants in dbsnp.
`bedtools intersect -u -a dbSnp153_chr1.bed.gz -b /users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.1.bed.gz`
Record the.. we got snp and delins. And there's rare/commonETC that has a definition.

The table we'll make is
1 - what percent of the catalog has any pVCF variant
2 - what percent of the catalog has any pVCF variant by size (snp, <5bp, 5, 50, 50+)

3 - what percent of the catalog has any dbsnp variant
	by types and by rare/common annotation.

If I intersect the catalog regions without pVCF or dbsnp common variants,
these are potentially regions which... something might be going on...
So I can run regioners... on these TRs and see if the...
... subset of non-variant TRs have more intersections to genes than the subset of variant TRs ...
And report in a supplementary table the set of TR regions which lack 'common' variants nearby genes
as well as information about their proximity to the genes


per-sample

hc_beds -> get variant counts
hc_beds -> subset catalog -> get variant counts
join/join


Going to stream down the two...
giant_coverage_allsamples.bed
and variants

bedtools intersect -u -f 1 -a catalog -b coverage | cut -f1-3 > samplename thing

Make a tsv of hc_bed, sample name

zip - 
I think the whole genome is covered...

cat all_hc_beds.txt | xargs bedtools multiinter -g ~/scratch/code/regione_rust/test_beds/grch38.genome.txt -names $(cat
all_names.txt) -i > giant_coverage_allsamples.bed


perm tests


~/code/regioners/target/release/regioners -g ~/code/regioners/test_beds/grch38.genome.txt --mask ../sf_trgeneprom/data/grch38.exclude_regions.bed -A TR_without_any_variant.bed -B ../sf_trgeneprom/data/grch38.protein_coding_transcript.bed -n 1000 --per-chrom --random circle -o TR_without_any_variant.regioner.json --threads 4

