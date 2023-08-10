For every site in the Tier1
    Collect the motif length from catalog (max)
    Collect the allele delta from the benchmark
    BoxPlot:
        X - binned motif length
        Y - allele delta distribution
    Repeat for observed allele deltas in pVCF non-HG002
    Do a hue of HG002 / Other samples
    Should show that the variety of variants in HG002 is
    representative of that shown over the population and therefore
    a good 'first pass' benchmark for tandem repeats

    Then we go into details about CODIS/Patho/MedRelevant Gene's 
    representation (most will be covered) for those that aren't
    explain why (limitations)

To make the new signed allele-delta, run
python new_allele_delta.py ~/code/adotto/benchmark/GIABTR_benchmark.6.26/GIABTR.HG002.benchmark.regions.bed.gz ~/code/adotto/benchmark/GIABTR_benchmark.6.26/GIABTR.HG002.benchmark.vcf.gz HG002 > n_ad.txt

This will recalculate output the max allele delta with contractions being negative and expansions being positive

Pct of TRs in SegDups remaining in benchmark:
bedtools intersect -a catalog.bed -b segdups.bed | wc -l
364589
bedtools intersect -a benchmark.bed -b segdups.bed | wc -l
127682
127682 of 364589 (35.0%)
