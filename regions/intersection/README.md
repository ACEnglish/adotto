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
- `counts_<output>` - input region entries annotated with number of variants and number of variant bases
- `filtered_<output>` - the counts file filtered to only regions containing non-SNP variants

And reports:
```
		v0.1		v0.3-dev
statistic       count   percent	count   percent
total regions   2232565 1	2170271 1
no variant      448124  0.2007	431781  0.1990
only a SNP      372144  0.1667	242294  0.1116
only SNPs       474209  0.2124	135780  0.0626
remaining       938088  0.4202	1360416 0.6268
```

Let's repeat this with the annotations we made previously
```
		v0.1		v0.3-dev
statistic       count   percent	count   percent
total regions   3298925 1	3503876 1
no variant      1600118 0.4850	1716435 0.4899
only a SNP      505514  0.1532	332505  0.0949
only SNPs       389598  0.1181	160201  0.0457
remaining       803695  0.2436	1294735 0.3695
```

And again with the unannotated regions
```
		v0.1		v0.3-dev
statistic       count   percent	count   percent
total regions   439538  1	428642  1
no variant      128123  0.2915	126221  0.2945
only a SNP      102119  0.2323	61672   0.1439
only SNPs       126488  0.2878	28007   0.0653
remaining       82808   0.1884	212742  0.4963
```

So it's interesting (promising) that our unannotated regions less frequently contain variants.

v0.3-dev ... We have a lot more regions 'remaining' in the unannotated. I gotta figure out what's happening here.

1. Adding these new regions (namely pbsv, trgt, and usc are expanding the boundaries. 
Collect these stats for the first slide... Actually hold off at this point.

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


Question 3
==========
Can we find expansions/contractions of the tr_annotations inside the variants?

The `filtered_variants_to_regions.txt` is now our new version of the tr_regions.bed. We'll use that to repeat the
'Defining Repeats' steps described in `../README.md` 
Then run TRF on the reference sequence of regions:

```bash
samtools faidx -r <(cat filtered_variants_to_regions.txt | awk '{print $1 ":" $2 "-" $3}') \
    ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa > tr_regions.fasta
trf409.linux64 tr_regions.fasta 3 7 7 80 5 5 500 -h -ngs > grch38.tandemrepeatfinder.txt
python ../scripts/trf_reformatter.py grch38.tandemrepeatfinder.txt final_something
bedtools sort -i final_something.bed | bgzip > final_something.bed.gz
tabix final_something.bed.gz
python ../scripts/tr_reganno_maker.py filtered_variants_to_regions.txt final_something.bed.gz > candidate_v0.3_anno.bed
```

Because we're going to be using the variants to filter these repeat annotations, we lower the min-score to 5 from 40
with the idea being we're more interested in sensitivity.

