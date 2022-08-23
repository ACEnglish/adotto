A Project-level VCF is created per-reference. Currently, only GRCh38 autosomes + X/Y are available.

# Versions:

## v0.1 - Initial version
(Click the badge to go to the download page)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6975244.svg)](https://doi.org/10.5281/zenodo.6975244)


### CHANGES:
* It exists

### Notes:
* Metadata for samples/individuals can be found [here](https://github.com/ACEnglish/adotto/tree/main/metadata)
* Assembly Sources:
  * HPRC
    * [Github](https://github.com/human-pangenomics/HPP_Year1_Assemblies)
  * Li
    * [FTP](ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/)
    * [Citation](https://scholar.google.com/scholar_lookup?author=S+Garg&author=A+Fungtammasan&author=A+Carroll&author=M+Chou&author=A+Schmitt&author=X+Zhou&title=Chromosome-scale%2C+haplotype-resolved+assembly+of+human+genomes&publication_year=2021&journal=Nat+Biotechnol&volume=39&pages=309-12)
  * Eichler
    * [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/)
    * [Citation](https://scholar.google.com/scholar_lookup?author=P+Ebert&author=PA+Audano&author=Q+Zhu&author=B+Rodriguez-Martin&author=D+Porubsky&author=MJ+Bonder&title=Haplotype-resolved+diverse+human+genomes+and+integrated+analysis+of+structural+variation&publication_year=2021&journal=Science&volume=372)
* Files provided per-chromosome as compressed binary VCF `.bcf.gz`
* To reduce file-size, FORMAT/FT field shortened to a single letter P[ass]|F[ail] and FORMAT/BPDP removed

### Paths:
Since the variants were uploaded per-chromosome, direct links to each of the files are provided below.

```
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr1.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr2.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr3.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr4.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr5.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr6.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr7.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr8.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr9.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr10.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr11.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr12.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr13.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr14.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr15.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr16.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr17.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr18.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr19.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr20.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr21.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chr22.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chrX.bcf.gz
https://zenodo.org/record/6975244/files/adotto_variants.grch38.sqoff.chrY.bcf.gz
```
