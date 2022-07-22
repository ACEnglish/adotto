

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
bedtools sort -i input.bed | bedtools merge | merged_bed_filter.py | bgzip > merged.bed.gz
tabix merged.bed.gz
```

Note that `merged_bed_filter.py` is the relative path to `scripts/merged_bed_filter.py`. 
This file removes regions that span fewer than 10bp or > 50kbp, or if they have an end position before than their start. It also creates a
`merging_stats.txt` with a simple summary statement.

Merging beds
============
Next, the bed files are merged

`stats` somewhere

Then I added a 75bp slop

Then I removed the variants >= .. 50kb

Then I removed variants within Nbp of reference gaps.

Then I ran TRF on the reference sequence of the remaining regions

Then translate that TRF output back to genomic coordinates and format to a bed/gz/tbi

Then I checked the input source beds against this set of regions to ensure that our new set
at least somewhat is representative of the input beds. For example, source ABC has a region
on chr:pos-end with motif GGG. Does this final produced bed have that same (or a similar) 
repeat description?

`stats` of the remaining regions



