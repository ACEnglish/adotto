
I want to clean/consoidate the giab raw calls
1) gotta take out the '#' header
2) add the file name to 4th column for something of an annotation
3) Then concat/merge


```bash
for i in GRCh38*.bed.gz
do
    name=${i%.bed.gz}
    bedtools sort -i <(zgrep -v "#" $i) | awk -v "name=$name" '{print $0 "\t" name}' | bgzip > fmt_$i
done
zcat fmt_GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz fmt_GRCh38_AllTandemRepeats_lt51bp_slop5.bed.gz \
     fmt_GRCh38_AllTandemRepeats_gt100bp_slop5.bed.gz fmt_GRCh38_AllTandemRepeats_gt10000bp_slop5.bed.gz \
     fmt_GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz \
     | bedtools sort | bgzip > giab_concat_input.bed.gz
```

Region without annotation
```
chr2    11274849        11274920
```

found in giab
```
chr2    11274874        11274895        GRCh38_AllTandemRepeats_lt51bp_slop5
```

grch38 sequence
```
>chr2:11274849-11274920
GTTACCTGAAACCATGGAGAGTACTGAATCCTATATATATATTATGTTTTTTCCTATACATACATATCTATG
```

Rerunning trf with score of 10

```
@seq
32 45 2 7.5 2 84 15 35 42 0 0 57 0.99 TA TATATATATATTAT GTTACCTGAAACCATGGAGAGTACTGAATCC GTTTTTTCCTATACATACATATCTATG
47 52 2 3.0 2 100 0 18 0 0 0 100 0.00 TT TTTTTT GTTACCTGAAACCATGGAGAGTACTGAATCCTATATATATATTATG CCTATACATACATATCTATG
56 66 4 2.8 4 100 0 33 54 18 0 27 1.44 ATAC ATACATACATA CTGAAACCATGGAGAGTACTGAATCCTATATATATATTATGTTTTTTCCT TCTATG
```

No variants in pVCF
