In total, the benchmark produces N variants ≥5bp and the replicate produces N (sup detail). The comparison metrics are
detailed in (sup detail). In summary, when running vcf_eval, we see an F1 (harmonic mean of precision and recall) of N,
Truvari bench F1 is N, and Truvari refine F1 of N. The reported performance metrics per-tool vary by variant length
(Figure 3). The performance of the benchmark according to vcf_eval drops off gradually after 5bp and fails to compare
variants over 1kbp in length. 10bp. Similarly, Truvari bench’s reported performance doesn’t begin until 50bp due to the
default size minimum parameter. However, Truvari refine reports higher performance across all size regimes and variant
types. 


1) Make a script to count the variants
```python count_variants.py ../../benchmark/GIABTR.HG002.benchmark.vcf.gz ../../benchmark/GIABTR.HG002.benchmark.regions.bed
142716
```

2) Generate the RTG stats (will need a script to subset to >=5bp in order to make the table comparison
3) Collect the alignment replicate bench/refine variant summaries
4) Make a script that will plot precision/recall by Type/Size for RTG, Bench,... Refine...
Refine I can use the laytr thing with phab


