Let's investigate the diversity of variants across samples.

Using [utmos](https://github.com/ACEnglish/utmos), we will count the number of unique variants per-sample.
We'll repeat this with just variants <10bp and those >=10bp and plot the percent of 'cdf' of variants over samples.

Commands
--------
```bash
utmos convert input.vcf.gz output.jl
utmos select -c 86 output.jl -o report.txt
```

A shortcut is to leverage utmos' "piece-wise" conversion by making per-chromosome jobs with `bash mk_convert_cmds.sh` to create per-chromosome conversion jobs 

Results
-------
![Utmos Plot](UtmosPlot.png)

We see that the variants >=10bp taper-off more slowly than the variants under 10bp. This suggests that there may be more diversity in larger alleles than smaller alleles.

### Limitations

#### samples vs individuals
The samples were not restricted to unique individuals, therefore the replicates may slightly confound these results. However, I believe this still serves as a reasonable estimation.

#### allele vs variant
This isn't actually measuring the allelic diverity but instead the diversity of variants. What I mean by this is that it is possible (likely) that larger alleles are have more variablility in their representations than those which are smaller. For example, an expansion/contraction of a tandem-repeat motif could be placed in multiple positions across the tandem-repeat sequence in the reference whereas a SNP generally has only a single position at which it can be placed. This could cause a single tandem-repeat allele to be appear non-shared between samples. Therefore, at least some of the 'increased diversity' we're observing in the larger alleles is due to poorly regularized/normalized/standardized variant representations.

