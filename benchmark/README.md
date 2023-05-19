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
3. Create the tiered/annotated bed file with `bash makeTiers.sh`

VCF creation
============

Step1 in tiering subset the HG002 HPRC assembly. Just use that for the VCF.

Extra Bed Annos
===============
- Add allele deltas using `scripts/allele_delta.py`
