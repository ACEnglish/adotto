# Random of an annotated region with the pacbio source
annotated entry
```
chr1    13601966        13602090        17      7.4     321     1.84    GGTGGCTCATGCCTTTG
```

regions part of region
```
chr1    13601936        13602107
```

pacbio merged bed etry
```
chr1    13601966        13602051
```

pacbio source entry
```
chr1    13601966        13602051        ID=chr1_13601966_13602051,STRUC=(GGTGGCTCATGCCTTTG)n
```

grch38 reference
```
>chr1:13601967-13602091
GGTGGCTCATGCCTTTGGGTGGCTCATGCCTTTGGGTGGCTCATGCCTTTGGGTGGCTCATGCCTTTGGGTGGCTCATGCCTTTGGGTGGCTCAGCACTTTGGGTGGCTCAGCACTTCGGGAGGC
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||***|||||||||||||****|||*|||*|||
GGTGGCTCATGCCTTTG                 GGTGGCTCATGCCTTTG                 GGTGGCTCATGCCTTTG                 GGTGGCTCATGCCTTTG
                 GGTGGCTCATGCCTTTG                 GGTGGCTCATGCCTTTG                 GGTGGCTCATGCCTTTG                 GGTGGC
									             <pacbio stops here (exact matches?)
```

# Missing 1
region without annotation
```
chr1    28636973        28637033
```
pacbio source entry
```
chr1    28636998        28637008        ID=chr1_28636998_28637008,STRUC=(TTTCA)n
```

grch38 sequence
```
>chr1:28636973-28637033
ACCCCATCTCTCAAAAAATACCAGTGTTTCATTTCACTCATCCTGCTGCTAACCCAACACA
                          |---||---|
```

pVCF
```
chr1    28636998        .       GTTTCA  G       60      PASS    AC=16;AN=26
```
So deletion of the TTCA is a fairly common allele

