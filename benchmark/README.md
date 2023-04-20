Creating a Benchmark
====================

After collecting TR regions and variants, we can create a benchmark of regions which are:

1. Confidently covered by an assembly
2. Confidently assessable by a pipeline

Inputs
======

- Adotto TRrv1.1 ([link](https://doi.org/10.5281/zenodo.7689784))
- Adotto pVCF ([link](https://doi.org/10.5281/zenodo.6975244))
- Dipcall HPRC HG002 coverage bed (provided by nist after manual curation)
- Dipcall HPRC HG002 VCF (provided by nist, `bcftools norm -m-any -N`)
- Adotto per-haplotype coverage beds
- GRCh38.chrX.PAR.bed pseudoautosomal regions of X
- GRCh38 genome bed file
- replicate_curated_states.txt which maps replicate states to tiers

```
inputs/
    adotto_TRregions_v1.1.bed.gz
    adotto_variants.grch38.sqoff.vcf.gz
    adotto_variants.grch38.sqoff.vcf.gz.tbi
    GRCh38_1kg_mainchrs.genome_bedtools.bed
    GRCh38.chrX.PAR.bed
    GRCh38_HG2-HPRC-20211005_dipcall-z2k.benchmark.bed
    HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz
    HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz.tbi
    HPRC.HG002.mat.cov.bed
    HPRC.HG002.pat.cov.bed
    replicate_curated_states.txt
```

Covered
=======

Create regions confidently covered by alignments of the HG002 HPRC assembly by running
```bash
scripts/covered_regions.sh
```

The main output file from this step is `covered/adotto_TRregions_v1.1_HPRC_HG002_covered.bed`, the subset of TRregions
which are confidently covered

Tiering
=======

Next, we need to use the pVCF HG002 replicates as well as the HPRC dipcall alignment replicate to begin splitting the
covered TRregions into Tier1 and Tier2 as well as adding some benchmarking annotations.

1. Subset HG002 and its assembly replicates from the pVCF. `bash hg002_subsets.sh`
2. Benchmark against the replicates `bash run_repl_bench.sh`
3. Run the makeTiers.sh pipeline (need to paste that together so its a single step. Actually, I need it to be 2 steps
   because there's the manual setting of replicate_states to the Tiers)
4. Tier1.bed, Tier2.bed are output into this directory (give better names)
5. Subset the pVCF to only HG002 (and maybe only the Tier1/Tier2 variants?)
6. Finished

Note: step 5 might be a place where I do the tr_identification_heuristics.py?
Though I'm not totally sure about that.

#.. I don't remember...
# I 
`adotto/strawman_chr20/clean_bench_regions/makeTiers.sh`
`adotto/strawman_chr20/clean_bench_regions/make_tiers.py`
`adotto/strawman_chr20/clean_bench_regions/split_tier1.py`

VCF creation
============

Subset to just HPRC HG002 and (optinally?) run `truvari anno trf`

`adotto/strawman_chr20/scripts/tr_identification_heuristics.py`

Intersect the assembly/alignment replicates and separate into tiers

Filter Tier1 regions without useful controls

Subset the pVCF to HG002 only.

Annotate with `truvari anno trf` over the well covered regions and run that through the script that sets the FILTER.


Should end up with:

GIABTR_HG002_Tier1.bed
GIABTR_HG002_Tier2.bed
GIABTR_HG002.vcf.gz
README.md
