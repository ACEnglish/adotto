In an attempt to clean up the TRregions for the HG002 benchmark, we're going to look for conflicting signals.

We have two other replicates of the HG002 assembly (Eichler/Li) as well as a separate alignment of the HPRC assembly.

We're going to benchmark over the Strawman bed regions by comparing the 'adotto' HPRC assembly with these three other
variant sets.

Then, we're going to compare the `refine.regions.txt` state reported by the different comparisons. Comparison of these
states will give us an idea of how good the region is for benchmarking. If all three agree, it's probably a good region.
If all three disagree, it could be a problematic region. There's many different state combinations across the three
comparisons. Each combination has been manually curated and assigned to a Tier1 or Tier2 signal. See states.txt.

After we have this state/tier assignment, we then run `make_tiers.py` to split our 'confidently covered' regions into
Tier1 and Tier2.

```bash
bash run_truvari_replicates.sh
bash run_truvari_alignments.sh
bash makeTiers.sh
```
