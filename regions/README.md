

We need to define what regions of the genome contain tandem-repeats (of interest).


Data setup
==========
We collected tandem repeat region bed files from various sources and merge them. See [DataDescription.md](DataDescription.md) for details

The sub directories of `data/` are each named for the contributor. Inside of those directories is `mk_merge.sh`
which will create the merged beds with something like

```bash
bedtools sort -i input.bed | bed_stats.py | bedtools merge | merged_bed_filter.py | bgzip > merged.bed.gz
tabix merged.bed.gz
```

Note that `bed_stats.py` and `merged_bed_filter.py` is the relative path to `scripts/merged_bed_filter.py`. 

`bed_stats.py` does a `pandas.Series.describe` on the input set of regions' span lengths and writes a summary to `input_spans.txt`

`merged_bed_filter.py` removes:
* regions not on chr1-22,X,Y
* regions that span fewer than 10bp or > 50kbp
* regions with an end position before than their start. 
It also creates a `merging_stats.json` with a simple summary statement.

Consolidating bed files
=======================
Next, the bed files are merged together to make the grand union using `data/mk_merge_grand.sh ref.genome`. This 
is almost the same as the component parts' `mk_merge.sh`. However, a ref.genome file needs to be provided (just 
`<chromName><TAB><chromSize>`)  This script creates `data/merged.bed.gz` as well as a `data/merged.slop25.bed.gz` 
The slop bed is a remerge of the `merged.bed.gz` after inflating all the regions by 50bp (25bp on each end).


Generating stats
================
All the intermediate stats files generated can be concatenated into a single tab-delimited table with
```
python scripts/consolidate_stats.py > data/region_stats.txt
```

Defining Repeats
================
At this point, we have broad regions where we expect to find tandem-repeats. Now we want to actually annotate the 
tandem repeat motifs, as well as their exact positions and their copy number in the reference.

* chrom - chromosome
* start - start position of the repeat
* end - end position of the repeat
* period - period size of the repeat
* copies - number of copes of the repeat in the reference
* score - alignment score
* entropy - entropy measure based on percent composition
* repeat - motif of the repeat

To do this, we'll start by running TRF on the regions. First, we extract the sequences from the regions

```bash
samtools faidx -r <(zcat tr_regions.bed.gz | awk '{print $1 ":" $2 "-" $3}') \
~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa > tr_regions.fasta
```

Then run TRF on the reference sequence of regions:
```bash
trf409.linux64 data/tr_regions.fasta 3 7 7 80 5 5 500 -h -ngs > data/grch38.tandemrepeatfinder.txt
```

TRF parsing
===========
Translate the trf output into a DataFrame and a bed file

```bash
python scripts/trf_reformatter.py data/grch38.tandemrepeatfinder.txt data/trf_annos
```
This creates `data/trf_annos.jl`, a DataFrame for analysis and `data/trf_annos.bed`.

TRF Intersection
================
Intersect the trf_annos back to the input sources for QC

```bash
python scripts/bed_intersection_stats.py
```
Note: run this inside the `regions/` directory. It has hard coded paths. Writes to `data/intersection.jl`

QC TRF
======
Now we want to put it all together.
We can QC the trf_annos.jl, and the `intersection.jl`, `region_stats.txt`, `trf_regions.bed`, etc

`notebooks/analysis.ipynb`

Unannotated regions
===================
Find `data/tr_regions.bed` without any `data/tr_annotations.bed` using
```bash
bedtools intersect -a tr_regions.bed.gz -b tr_annotated.bed -c | awk '$4 == 0' > unannotated_regions.bed
```

Inside `manual_inspection/` are notes of looking at annotated and unannotated regions when intersected back to the
sources.
Tips - commands for pulling stuff:
```
tabix some.bed.gz chr:start-end # Query bedfiles
samtools faidx reference.fa chr:start-end # Query for sequence
bcftools view -r chr:start-end some.vcf.gz # Query variants
```
Note `samtools faidx` is not the same coordinate system as `tabix`. But the bed files are 0-based and vcfs are 1-based,
both half-open(?).


Creating the intermediate annotations
=====================================

As of right now, we have tr regions (simple chr:start-end) and annotations (coordinates + other info). Since there can
be multiple annotations per-region, it makes sense to combine these two files into one to facilitate analysis. We'll
design this file for use with `truvari anno trf`, but it can be adopted to other uses. The final file format will be a
tab-delimited file with columns:

* chrom - tr_region chromosome
* start - tr_region 0-based start position
* end - tr_region 0-based end position
* annos - json of annotations

The json of annotations will be a list of dicts with key:values listed above in "Defining Repeats".

To create this file, simply run:

```bash
python scripts/tr_reganno_maker.py intersection/data/tr_regions.bed intersection/data/tr_annotations.bed.gz \
    | bedtools sort -i - | bgzip > candidate_v1.bed.gz
tabix candidate_v1.bed.gz
```

Creating RepeatMasker annotations
=================================

These regions may have 'contamination' from non-TR repeats which TRF will discover (e.g. tandem ALUs). Therefore, we
next want to analyze the regions' sequences for non-TR repeats so that we can potentially filter them. We accomplish
this by

Run `scripts/fasta_splitter.py` on tr_regions.fasta, then RepeatMasker the parts
```bash
mkdir -p jobs
for i in part*fasta
do
    
    echo "#!/bin/bash" > jobs/${i}.sh
    echo RepeatMasker -pa 16 -qq -e hmmer -species human \
        -lcambig -nocut -div 50 -no_id -s $i >> jobs/${i}.sh
done
```
Then make a simple to parse dataframe with `scripts/repmask_to_joblib.py`


Improving the annotation
========================

The current state of the intermediate TR regions lacks other useful information.
We add them in by running 

```bash
python scripts/annotation_improver.py \
	data/candidate_v1.0.bed \
	grch38.fasta \
	repmask.jl \
	pathorepeats.bed \
	codisrepeats.bed
	> adotto_TRregions_v1.0.bed
```

Then we use `scripts/make_gene_lookup.py` and `scripts/append_gene.py` to add the gene columns
