regions+variants intersection
=============================

Now that we have the candidate regions/annotations and a large pVCF of variants, we can answer some questions:

1. How many candidate regions have no observed variation?
2. Of the candidate regions with variation, what percent of the variants by count and bases effected are contained
   within?
3. Can we find expansions/contractions of the tr_annotations inside the variants?


Question 1:
===========

Run the following:
```bash
python variant_region_intersection.py tr_regions.bed.gz adotto_variants.grch38.sqoff.bcf.gz variants_to_regions.txt
```
Note, to prevent large SVs which span regions from altering our counts, the variants' start/end boundaries must be
within the region's boundaries.

This creates two files:
- `counts_variants_to_regions.txt` - input regions.bed entry annotated with number of variants and number of variant bases
- `filtered_variants_to_regions.txt` - the counts file filtered to only regions containing any non-SNP variants

And reports:
```
statistic       count   percent
total regions   2232565 1
no variant      448124  0.2007
only a SNP      372144  0.1667
only SNPs       474209  0.2124
remaining       938088  0.4202
```

Let's repeat this with the annotations we made previously
```
statistic       count   percent
total regions   3298925 1
no variant      1600118 0.4850
only a SNP      505514  0.1532
only SNPs       389598  0.1181
remaining       803695  0.2436
```

And again with the unannotated regions
```
statistic       count   percent
total regions   439538  1
no variant      128123  0.2915
only a SNP      102119  0.2323
only SNPs       126488  0.2878
remaining       82808   0.1884
```

So it's interesting (promising) that our unannotated regions less frequently contain variants.

Question 2:
===========
Of the candidate regions with variation, what percent of the variants by count and bases effected are contained
within?

Let's start by getting the 'negation' of our tandem repeat regions
```bash
bedtools genomecov -i counts_variants_to_regions.txt -g genome.fa.fai -bga \
	| awk '$4 == 0' > uncovered_v0.1_regions.bed
```
And repeat for `filtered_variants_to_regions.txt`, which is probably going to be our `tr_regions_v0.2.bed`

Then we can do some pandas math to see:

```
version	region	var	count		bases
v0.1	TR	all	24922545	494028725
v0.1	non-TR	all	99779718	584875735
v0.1	TR	sv	556398		429845954
v0.1	non-TR	sv	186699		477883792
v0.2	TR	all	21882069	490988249
v0.2	non-TR	all	102828813	594648991
v0.2	TR	sv	556398		429845954
v0.2	non-TR	sv	190319		484583436

For context:
v0.1 TR genome coverage:	238052458bp	7.44%
v0.2 TR genome coverage:	121788538bp	3.81%
```
So, 3.81% of the genome which the v0.2 TR regions cover contains
- 17.5% of all variants by count
- 45.2% of all variants by bases effected
- 75.5% of SVs by count
- 47.0% of SVs by bases effected
