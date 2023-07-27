GIAB TR Benchmark HG002 GRCh38
==============================

See below for information on how to use this benchmark.

To see descriptions of the steps used to create this benchmark, see [the github](https://github.com/ACEnglish/adotto/tree/main/benchmark)

Changes
=======
6.26
- Altering of Tier1/Tier2 definitions
- Adding Maternal/Paternal haplotype allele delta columns to benchmark regions bed
- Adding `laytr giabTR` report supplementary files and documentation
- README updated with Truvari v4.1-dev usage

4.20
- Correction of 'confidently covered' chrX,chrY with respect to pseudoautosomal regions
- Creation of Tiers [slides](https://github.com/ACEnglish/adotto/blob/main/slides/GIABTR_English_April182023.pdf)
- Simplified files

Usage
=====
To benchmark against the regions, use the development branch of Truvari >= 4.1-dev

```
git clone https://github.com/ACEnglish/truvari.git
cd truvari
python3 -m pip install .
```

The two commands you'll need are:
```
truvari bench -b GIABTR.HG002.benchmark.vcf.gz \
	-c <your_calls.vcf.gz> \
	--includebed GIABTR.HG002.benchmark.regions.bed.gz \
	--sizemin 5 --pick ac -o bench_result/
truvari refine --use-original-vcfs --reference grch38.fasta bench_result/
```

If you want to check only regions analyzed by your caller (i.e. not penalizing your caller for not analyzing all TR regions), use the command:
```
truvari refine --regions <your_regions.bed> --use-original-vcfs --reference grch38.fasta bench_result/
```

See [truvari wiki](https://github.com/ACEnglish/truvari/wiki) for details on the outputs created and information on
other parameters available for `bench` and `refine`

For detailed stratifications, install [laytr](https://github.com/ACEnglish/laytr) via:

```
git clone https://github.com/ACEnglish/laytr.git
cd laytr
python3 -m pip install .
```

You can then generate an html report with performance for different stratifications via:
```
laytr giabTR --regionsummary bench_result/refine.regions.txt\
	--includebed GIABTR_benchmark.6.26/GIABTR.HG002.benchmark.regions.bed.gz \
	--som GIABTR_benchmark.6.26/adotto_TRv1.1_4mers.som \
	--somap GIABTR_benchmark.6.26/adotto_TRv1.1_4mers.map
	--trcatalog adotto_TRregions_v1.1.bed.gz
	--output your_report.html
```

Where `bench_result/refine.regions.txt` is the per-region information from the `truvari bench/refine` steps performed
above and `adotto_TRregions_v1.1.bed.gz` is the TR catalog found [here](https://github.com/ACEnglish/adotto/blob/main/regions/DataDescription.md).

Files
=====
GIABTR.HG002.benchmark.regions.bed.gz
-------------------------------------
The tandem repeat benchmark. Each region is a TR from adotto TR regions v1.1 which is confidently covered by the HPRC
HG002 assembly
Columns:
- chrom
- start
- end
- Tier - 1 for confident variant state and confidently assessable. 2 for everything else
- Replicates - underscore-delimited triple of region benchmarking states when compared to two HG002 assembly technical replicates and an alignment replicate
- Variant State flag - bitflag of variants found within region
-- 0x1 : HG002 >=5bp variant
-- 0x2 : HG002 <5bp variant
-- 0x4 : Other >=5bp variant
-- 0x8 : Other <5bp variant
- Spanned sequence entropy
- Maternal haplotype variant length sum
- Paternal haplotype variant length sum

GIABTR.HG002.benchmark.vcf.gz
-----------------------------
The HG002 only subset of [adotto v0.1 variants](https://github.com/ACEnglish/adotto/blob/main/variants/DataDescription.md)


adotto_TRv1.1_4mers.[map|som]
-----------------------------
Supplementary files for use with `laytr giabTR` report creation ([details](https://github.com/ACEnglish/laytr))
