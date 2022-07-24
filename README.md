Project Adotto
==============

This repo holds the analysis scripts/notes/summaries for the Genome in a Bottle Tandem Repeat Benchmark.


For now, we'll be assuming all analysis is using GRCh38 as the reference

Structure
=========
There are a few main sub-parts to this project. Each is contained in a sub-directory. Raw data that's too large to be
kept on github will be made available and  documented such that a user can find it and know where to place it within
a clone of this repo in order to run sub-parts of the analysis.

* regions - Identification of Tandem-Repeat regions of a reference
* variant_calling - Calling variants from long-read haplotype resolved assemblies
* realignment - Joint-realignment of tandem-repeat regions of a reference
* benchmarking - Scripts and results from benchmarking work
* manuscript - Summary and plotting workflows for the publication
* metadata - Data descriptor files (e.g. download paths of inputs used or sample ancestry information)
* slides - GIAB team meeting slides

Quick Start
===========

The first few iterations of this project will not focus on software usability. Instead, we're only trying to collect
verbose documentation. As the project matures, we'll continually improve code so that the work can more easily be
recreated and/or improved upon. 

Please feel free to open Issues in this repository for any questions - software related or otherwise.

Data
====
Data can be packaged up and put in a tar with the same directory structure using:
```bash
find base -name "*.txt" | tar czvf base.tgz -T -
```
Then, unzipping the data from the same directory will put the files in-place

Progress
========
Rough estimates of the progress for each analysis section.  Last updated  `Fri Jul 22 02:03:49 CDT 2022`
| Step            |  %  |
| --------------- | --- |
| regions         | 80% |
| variant_calling | 95% |
| realignment     | 20% | 
| bench creation  | 5%  |
| benchmarking    | 0%  |
