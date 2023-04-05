In an attempt to clean up the TRregions for the HG002 benchmark, we're going to look for conflicting signals.

We have two other replicates of the HG002 assembly (Eichler/Li) as well as a separate alignment of the HPRC assembly.

We're going to benchmark over the Strawman bed regions by comparing the 'adotto' HPRC assembly with these three other
variant sets.

If we find sites where Eichler/Li have False Positives, that could be indicative of collapsed hets, variants
missed by the HPRC assembly (unlikely), or alignment ambiguities pulling Eichler/Li variants into the regions.

If we find FP/FN when comparing to the 'dipcall' HPRC assembly, that's evidence of alignment problems or Truvari's
inability to make an accurate comparison.


```bash
bash makeTiers.sh
```
