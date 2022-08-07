Benchmarking QC
===============

We now want to check the quality of the variants we called by comparing them to truth-sets. We have two
truth-sets on GRCh38 from GIAB:

- [CMRG](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz)
- [CMRG bed](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed)
- [CMRG smallvar](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz)
- [CMRG smallvar bed](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed)
- [TrioHifiAsm](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/HPRC-HG002.cur.20211005/HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz)
- [TrioHifiAsm bed](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/HPRC-HG002.cur.20211005/HPRC-cur.20211005-align2-GRCh38.dip.bed)

Note that the GIAB HPRC Dipcall uses the same long-read haplotype resolved assembly as the pVCF's HG002 sample.

We'll run rtg and truvari over our 3 replicates of (HG002) NA24385 against CMRG and the TrioHifiAsm. We'll then make notebooks
to analyze and report the results.


Setting up inputs
-----------------
Using rtg version 3.12.1, we format the index via:
```bash
rtg format -o GRCh38_1kg_mainchrs.sdf GRCh38_1kg_mainchrs.fa
```

While truvari has functionality to pick a single sample from a VCF and to ignore uncalled (reference-homozygous) entries from a
sample, rtg does not appear to have the same. Therefore we need to also pull out the comparison samples from our pVCF:
```bash
for i in HG002 NA24385 li:NA24385
do 
  bcftools view -s $i GRCh38.variants.sqoff.vcf.gz | vcf-subset -e | bgzip > ${i}_comp.vcf.gz
  tabix ${i}_comp.vcf.gz
done
```

Finally, truvari expects single-allele vcfs whereas the THFA holds multi-allelics. Therefore, we must split those
variants:
```bash
bcftools norm -m-any -N HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz \
	-O z -o HPRC-cur.20211005-align2-GRCh38.dip.singlealleles.vcf.gz
```

Running benchmarking
--------------------
Edit the hard coded paths at the top of `scripts/run_bench.sh` to point to the appropriate files and then run the
scripts to create all the `results/`.


Consistency
-----------
Manually inspecting the `results/truvari_thfa*/summary.txt` shows that there's lower recall than I would expect with
~10k FN per-replicate in the thfa 
```bash
truvari consistency -j results/truvari_thfa_*/fn.vcf > results/missing_stats.json
```
The report shows that of the 11,480 FNs across replicates, 9,912 are missed by all three replicates. This suggests
there's a problem with the alignments. Let's pull out the shared FNs and query the coverage annotation.

No Coverage Regions
-------------------
Concatenate, sort, and calculate how many haplotypes have no aligned coverage.
```bash
awk '{print $2 "\n" $3}' coverage_files.txt \
	| xargs cat \
	| awk '$4 == 0' \
	| bedtools sort -i - \
	| bedtools genomecov -i - -g genome.bed -bga | awk '$4 == 172' > never_covered.bed
```
Let's see how many FNs are contained in these regions.. 698 FNs are in the `never_covered.bed`. Let's try this again
with only the HG002 assemblies.. And it's only up to ~740.

Perhaps not using multimatch is the problem?
...
```
for i in truvari_thfa_*/fn.vcf; do bcftools view -i "Multi != 1" $i | grep -v "#" | cut -f1-5; done | sort | uniq -c | sort -nr | awk '{$1=$1};1' | cut -f1 -d\  | sort | uniq -c
    886 1
    865 2
   4720 3
for i in truvari_thfa_*/fn.vcf; do bcftools view -i "Multi == 1" $i | grep -v "#" | cut -f1-5; done | sort | uniq -c | sort -nr | awk '{$1=$1};1' | cut -f1 -d\  | sort | uniq -c
    822 1
    431 2
   4525 3
```
So it's 50/50 if the always missed variants are FNs because of multimatching.

Okay, so I need to Check if they're multicovered or whatever. So I'll use the pVCF anno altered so that I can 
.. But multi-covered isn't
