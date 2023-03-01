# Versions:

## v1.1
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7689784.svg)](https://doi.org/10.5281/zenodo.7689784)

- Updated pathogenic repeats. 54 pathogenic regions are unchanged, 2 have been changed, and 6 added.
```
changed: (is_now -> was)
	NOTCH2NLA -> NOTCH2NLC
	NOTCH2NLC -> NOTCH2NL
added:
	EIF4A3, PRNP, TBX1, PRDM12, DMD, ZIC3
```
- Renamed annotations' "motif" key back to "repeat" for `truvari anno trf` compatibility.
- hom_span column normalized as hom_pct - percent of bases in the regions annotated as homopolymers

## v1.0
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521434.svg)](https://doi.org/10.5281/zenodo.7521434)

### CHANGES: 
* Removed homopolymer annotations
* Simplified overlapping annotations
* Added new columns that describe properties of the region 

### File Structure:
| Column        | Definition                                                                                     |
|---------------|------------------------------------------------------------------------------------------------|
| chr           | Chromosome of the region                                                                       |
| start         | Start position of the region                                                                   |
| end           | End position of the region                                                                     |
| ovl_flag      | overlap categories of annotations inside the region                                            |
| up_buff       | number of bases upstream of the first annotation's start that are non-TR sequence              |
| dn_buff       | number of bases downstream of the last annotation's end that are non-TR sequence               |
| hom_span      | number of bases of the region found to be homopolymer repeats                                  |
| n_filtered    | number of annotations removed from the region                                                  |
| n_annos       | number of annotations remaining in the region                                                  |
| n_subregions  | number of subregions in the region                                                             |
| mu_purity     | average purity of annotations in region                                                        |
| pct_annotated | percent of the region's range (minus buffer) annotated                                         |
| interspersed  | name of interspersed repeat class found within region by RepeatMasker v4.1.4                   |
| patho         | name of gene affected by a pathogenic tandem repeat in region                                  |
| codis         | name of CODIS site contained in region                                                         |
| gene_flag     | gene features intersecting region (Enseml v105)                                                |
| biotype       | comma separated gene biotypes intersecting region (Enseml v105)                                |
| annos         | JSON of TRF annotations in the region (list of dicts with keys: motif, entropy, ovl_flag, etc) |

The `annos` JSON is simple key:values of:
* chrom - chromosome
* start - start position of the repeat
* end - end position of the repeat
* period - period size of the repeat
* copies - number of copes of the repeat in the reference
* score - alignment score
* entropy - entropy measure based on percent composition
* motif - motif sequence of the repeat
* purity - Sequence similarity of  `motif*copies` against annotation’s reference span


## v0.3 - More Regions
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7226352.svg)](https://doi.org/10.5281/zenodo.7226352)

### CHANGES:
* Added new annotations sources from:
  * [TRGT](https://github.com/PacificBiosciences/trgt/tree/main/repeats) - Both full regions and pathogenic
  * [pbsv](https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed)
  * [Vamos](https://zenodo.org/record/7155334/)
* See [slides](https://github.com/ACEnglish/adotto/blob/main/slides/GIABTR_English_October172022.pdf) for details
* Same file structure as v0.2


## v0.2 - Useable version
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7013709.svg)](https://doi.org/10.5281/zenodo.7013709)

Much simpler and much smaller data

### CHANGES:
* Using a pVCF created from 172 haplotype-resolved long-read assemblies, we removed any TR region which had no observed
  non-SNP variant. This removed 58% of the regions.
* With fewer regions, we simplified the data's structure and now store everything in a single tab-delimited bed-like
  file with columns:

### File Structure:
Columns:
* chrom - chromosome of TR region
* start - 0-based start position of TR region
* end - 0-based end position of TR region
* annos - Json containing a list of TRF annotated repeats with structure:
  * chrom - chromosome
  * start - start position of the repeat
  * end - end position of the repeat
  * period - period size of the repeat
  * copies - number of copes of the repeat in the reference
  * score - alignment score
  * entropy - entropy measure based on percent composition
  * repeat - motif of the repeat

Example:
```
chr11   11605859        11605958        [{"chrom": "chr11", "start": 11605912, "end": 11605927, "period": 2.0, "copies": 8.5, "score": 41, "entropy": 0.99, "repeat": "CA"}, {"chrom": "chr11", "start": 11605930, "end": 11605941, "period": 6.0, "copies": 2.0, "score": 36, "entropy": 1.92, "repeat": "AGCTTC"}]
```
### Notes
The single file of regions paired with their annotations allows much easier parsing/usage.
A custom parser can easily be built by splitting each line on tabs and using a json parser on the 4th column.
A parser has already been built into [truvari](https://github.com/ACEnglish/truvari) and can be used via:

```python
from truvari.annotations.trf import iter_tr_regions
# generator for every region
adotto_all_regions = iter_tr_regions("adotto_TRannotations_v0.2.bed.gz")
# fetch regions in a region
adotto_fetch_regions = iter_tr_regions("adotto_TRannotations_v0.2.bed.gz", region=("chr17", 10350000, 10360000))
```

Additionally, the file can be queried with [tabix](http://www.htslib.org/doc/tabix.html)
```bash
tabix adotto_TRannotations_v0.2.bed.gz chr11:11600000-11606000
```

## v0.1 - Initial version
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6930202.svg)](https://doi.org/10.5281/zenodo.6930202)

Inside of data_list.txt is the relative paths to data generated by this step.

We bundle these files using tar such that they can be recovered in-place into a copy of the code for reanalysis and to
keep things organized.

The tarball can be created with the command:

```bash
cat data_list.txt | tar czvf adotto_regions_data_<version>.tgz -T -
```

The tarball can be placed inside this directory and extracted in-place via:

```bash
tar xzvf adotto_regions_data_<version>.tgz
```
Note! the data_list.txt is put into .gitignore to help keep `git status` clean 

### CHANGES:
* It exists

### Notes:
List of files and their descriptions:
* `data/tr_regions.bed.gz` - Final set of tandem-repeat regions for analysis
* `data/tr_annotated.bed.gz` - TandemRepeatFinder annotations over the tandem-repeat regions
* `data/unannotated_regions.bed.gz` - tr_regions.bed.gz which have no accompanying tr_annotated.bed.gz entries
* `data/merged.slop25.bed.gz` - Merged calls from the sources with 25bp of slop added to each end
* `data/giab/giab_concat_input.bed.gz` - Raw input provided by GIAB from [ftp](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/)
* `data/giab/merged.bed.gz` - GIAB regions merged
* `data/baylor/grch38.simpleRepeat.truvari.bed.gz` - Raw input provided by baylor from [UCSC Simple repeat annotations](https://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=varRep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema)
* `data/baylor/merged.bed.gz` - Baylor regions merged
* `data/pacbio/repeat_catalog.hg38.bed.gz` - Raw input provided by pacbio from [Illumina](https://github.com/illumina/Repeatcatalogs)
* `data/pacbio/merged.bed.gz` - Pacbio regions merged
* `data/ucsd1/ensembleTR_loci_list.bed.gz` - Raw input provided by UCSD from ensemble
* `data/ucsd1/merged.bed.gz` - UCSD1 regions merged
* `data/ucsd2/GIAB_adVNTR_short_VNTR_regions.bed.gz` - Raw input provided by UCSD from ???
* `data/ucsd2/merged.bed.gz` - UCSD2 regions merged
