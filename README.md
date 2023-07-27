Project Adotto
==============

This repo holds the analysis scripts/notes/summaries for the Genome in a Bottle Tandem Repeat Benchmark.

The Benchmark
=============

## 💾 Download
The current beta-release (v6.26) of the benchmark can be found
[here](https://drive.google.com/file/u/1/d/1fbOIUr7c4wdJ_Ede7hJlhAq-keX-yyuO/view?usp=sharing).

## 📜 README
Included in the above tar-ball is a `README.md` which can be seen
[here](https://github.com/ACEnglish/adotto/tree/main/benchmark/benchmark/GIABTR_benchmark.6.26/README.md).
This README has details on the files contained in the benchmark as well as instructions on how to compare your
caller using Truvari. Note that Truvari is currently under active development as a part of the GIAB TR project 
so a manual installation of the latest development branch is recommended. The README explains this in detail.


Other Data
==========
The full adotto tandem-repeat catalog and pVCF of 86 haplotype-resolved assemblies created as part of this work
are available via [zenodo.org](zenodo.org). Details of the files as well as download links are below.

| Dataset | Zenodo | Current Version | Details |
| ------- | ------ | --------------- | ------- |
| Regions | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7689784.svg)](https://doi.org/10.5281/zenodo.7689784) | v1.1 | [Link](regions/DataDescription.md) | 
| Variants | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6975244.svg)](https://doi.org/10.5281/zenodo.6975244) | v0.1 | [Link](variants/DataDescription.md) |

This repository
===============
There are a few main sub-parts to this project. Each is contained in a sub-directory. Raw data that's too large to be
kept on github will be made available and documented such that a user can find it and know where to place it within
a clone of this repo in order to run sub-parts of the analysis.

* slides - GIAB team meeting slides
* manuscript - Summary and plotting workflows for the publication
* metadata - Data descriptor files (e.g. download paths of inputs used or sample ancestry information)
* regions - Identification of Tandem-Repeat regions of a reference
* variants - Calling variants from long-read haplotype resolved assemblies
* benchmark - Scripts for consolidating regions and variants to create the benchmark
