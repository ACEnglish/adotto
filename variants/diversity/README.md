Let's investigate the diversity of variants across samples.

Using [utmos](https://github.com/ACEnglish/utmos), we will count the number of unique variants per-sample.
We'll repeat this with just variants <10bp and those >=10bp and plot the percent of 'cdf' of variants over samples.

The hypothesis is that variants greater than 10bp will 'taper-off' more slowly than those less 10bp, which might suggest
that there is more 'diversity' between samples in larger events.

Commands
--------
```bash
utmos convert input.vcf.gz output.jl
utmos select -c 86 output.jl -o report.txt
```

A shortcut is to leverage utmos' "piece-wise" conversion by making per-chromosome jobs with `bash mk_convert_cmds.sh` to create per-chromosome conversion jobs 

Results
-------
See `UtmosAnalysis.ipynb`
