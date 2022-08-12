pVCF
====

A pipeline for turning long-read assemblies into project level VCFs.

# Step 0 - requirements

Below is a list of tools/libraries used by the code. Setup your environment to have them all available 

!!TODO - update all these link

Requirements:
* minimap2 - v2.24 [link](https://github.com/lh3/minimap2)
* paftoools - packed with minimap2
* bcftools - v1.12
* vcftools - 0.1.16
* bedtools - v2.30.0 [link](https://github.com/arq5x/bedtools2)
* bgzip/tabix
* python3
  * truvari - v3.4+ [link](https://github.com/ACEnglish/truvari)

# Step 1 - organize the input assembly files

Inside of the `assemblies` folder, there's a file named `relative_path_metadata.txt`.
This holds the metadata and relative file paths to the input long-read haplotype resolved assemblies.
The columns are:
* sample - name of the sample
* project - name of the project which produced the assembly
* haplotag - identifier of which haplotype the fasta file contains
* rel_path - relative path of the assembly

This file helps with batch job submission.

For HPRC assemblies (in-progress), we downloaded the AGC file and extracted per-sample as needed


# Step 2 - map each haplotype

Each assembly is a single haplotype of a sample. `scripts/map_haplo.sh` will align a fasta to a reference, call
variants, create coverage beds, generate basic qc stats, and produce an annotated VCF. This summary is explored in `notebooks/HaplotypeCoverage.ipynb`

usage
-----
```bash
bash scripts/map_haplo.sh input.fasta reference.fasta sample_name out_dir
```
Additionally, there are hard-coded parameters towards the top of the script that can be changed.
* params - parameters provided to minimap
* threads - number of threads to use
* min_qual - minimum quality score of alignments used for variant calling or coverage bed creation

Note: the other downstream steps expect `out_dir` to be a sub-directory of `initial_alignments` e.g.
`initial_alignment/giab_HG002_mat`

Note: the reference.fasta.fai needs to be present

outputs
-------
* aln.paf - the minimap2 raw alignments
* aln.paf.stats.txt  - paftools stats on the aln.paf
* aln.bed - converted aln.paf to bed
* cov.bed - bedtools genomecov of the aln.bed
* single_cov.bed - subset of cov.bed of only singly covered genomic regions
* aln.vcf.gz[.tbi] - paftools called variants
* aln.covanno.vcf.gz[.tbi] - aln.vcf.gz annotated with cov.bed with FORMAT/BPDP and FILTER==COV for any variants not
  singly covered.

batching
--------
Use `scripts/mk_initial_alignment_jobs.sh` to create all the mapping commands at once. This parses the
`assemblies/relative_path_metadata.txt` to make one `job/aln_${project}_${sample}_${haplotag}.sh` script
per-haplotype. These scripts can then be easily submitted to a cluster.

Note: the reference is hard-coded into this script.

# Step 2.5 - create mapping stats

One or more directories produced in step 2 can be parsed with `scripts/hap_mapstats.py` and turned into a dataframe saved as
`hap_mapstats.jl`. The directory's name is used as the index for the results.

# Step 3 - combining haplotypes per-sample

Once a haplotype-pair is through initial alignment, two VCFs need to be combined to make the samples' VCF. Use
`scripts/merge_haps.sh` to do the combining. This step will use `bcftools merge` to combine the two VCFs, matching
exactly identical alleles, use the `scripts/single_sample_cov.py` to consolidate the SAMPLE information, add coverage
relevant FORMAT fields, add INFO/SVLEN,INFO/SVTYPE annotations, and remove the minimap INFO/Q.\* fields

usage
-----
This uses
```bash
bash scripts/merge_haps.sh hap1.vcf.gz hap1.bed hap2.vcf.gz hap2.bed sample output.vcf.gz
```

* hap\*.vcf.gz - the aln.covanno.vcf.gz from Step 1
* hap\*.bed - the cov.bed from Step 1
* sample - name of the sample to be put into the VCF
* output.vcf.gz - destination of the output VCF file

I organized the outputs into e.g. hapo_merged/eichler_HG00864.vcf.gz

# Step 4 - merging between samples

First, we make the header
```bash
bcftools merge -m none hapo_merged/*.vcf.gz -o pVCFs/GRCh38.variants.header.vcf --print-header --force-samples
```
Our run had redundant sample names which were manually altered inside of the output header e.g. NA24385 is assembled by
both eichler and li, therefore the li was renamed in the new header li:NA24385.

Next, we merge the vcfs
```bash
bcftools merge -m none hapo_merged/*.vcf.gz --use-header pVCFs/GRCh38.variants.header.vcf \
	| bcftools view -S pVCFs/sample_order.txt -o pVCFs/GRCh38.variants.vcf.gz -O z
```
The sample_order.txt is just the list of the sample names so that we can control their order in the pVCF.

# Step 5 (optional) - create shards
At this point, the GRCh38.variants.vcf.gz is pretty big. To speed up the squaring-off step, we divide it into shards

```bash
truvari divide -m 1000000 pVCFs/GRCh38.variants.vcf.gz pVCFs/parts/
```

# Step 6 - squaring off
Naturally, there are a lot of missing (`./.`) genotypes inside the pVCF. We can square it off with the tool

```bash
python scripts/generate_coverage_tree.py ../metadata/coverage_files.txt
```
The coverage_files.txt just organizes the pairs of cov.bed files per-sample
This script creates three files that make it easier to look up coverage across all the samples for squaring-off the VCF.
* annotree.jl - joblib saved intervaltrees per-chromosome with keys that link to the coverage entries
* coverage.jl - joblib saved pandas dataframe with keys from the intervaltrees and actual coverage depths

```bash
python scripts/annotate_pvcf_cov.py GRCh38.variants.vcf.gz annotree.jl coverage.jl | bgzip > GRCh38.variants.sqoff.vcf.gz
```

Note: This is slow. I ran one job for each `pVCFs/parts/*.vcf.gz` and they each took approximately an hour.

# Step 7 (optional) - combine the parts
If Step 4 was performed, and you used truvari v3.4-dev+ such that the divide output VCF names have `ls`-sortable vcfs,
recombine the parts using
```bash
vcf-concat pVCFs/sqoff_parts/*.vcf.gz | bgzip > pVCFs/GRCh38.variants.squareoff.vcf.gz
```

Otherwise, you'll have to do `vcf-concat | vcf-sort | bgzip > output.vcf.gz`
