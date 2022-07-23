

We need to define what regions of the genome contain tandem-repeats (of interest).


Data setup
==========
We collected tandem repeat region bed files from various sources and merge them

* [GIAB annotations](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/)
* [UCSC Simple repeat annotations](https://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=varRep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema)
* Pacbio - To be specified
* UCSD1 - To be specified
* UCSD2 - To be specified

The sub directories of `data/` are each named for the contributing source. Inside of those directories is `mk_merge.sh`
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
Next, the bed files are merged together to make the grand union using `data/mk_merge_grand.sh ref.genome`. This is almost the same as
the component parts' `mk_merge.sh`. However, a ref.genome file needs to be provided (just `<chromName><TAB><chromSize>`) 
This script creates `data/merged.bed.gz` as well as a `data/merged.slop25.bed.gz` The slop bed is a remerge of the
`merged.bed.gz` after inflating all the regions by 50bp (25bp on each end).

Generating stats
================
All the intermediate stats files generated can be concatenated into a single tab-delimited table with
```
python scripts/consolidate_stats.py > data/region_stats.txt
```

Reference Gaps
==============
Remove regions within 5kbp of reference gaps obtained from the 'Mapping and Sequencing>Gap' track from 
[UCSC](https://genome.ucsc.edu/cgi-bin/hgTables)

```bash
bash ../scripts/remove_gaps.sh HumanGRCh38.mapping.gap.bed.gz grch38.genome merged.slop25.bed.gz
```
This removes 1769 regions and puts the remaining ones in `data/tr_regions.bed.gz`

TRF on regions
==============
We now want to (re)create all the tandem repeat annotations. First, we extract the sequences from the regions

```bash
samtools faidx -r <(zcat tr_regions.bed.gz | awk '{print $1 ":" $2 "-" $3}')
~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa > tr_regions.fasta
```

Then I ran TRF on the reference sequence of regions:
```bash
trf409.linux64 data/tr_regions.fasta 3 7 7 80 5 40 500 -h -ngs > data/grch38.tandemrepeatfinder.txt
```

TRF to DF
=========
Translate the trf output into a DataFrame with

```bash
python scripts/trf_reformatter.py data/grch38.tandemrepeatfinder.txt data/trf_annos_df.jl
```
This creates `data/trf_annos_df.jl`

QC TRF
======
!! bookmark


Then I checked the input source beds against this set of regions to ensure that our new set
at least somewhat is representative of the input beds. For example, source ABC has a region
on chr:pos-end with motif GGG. Does this final produced bed have that same (or a similar) 
repeat description?

`stats` of the remaining regions



