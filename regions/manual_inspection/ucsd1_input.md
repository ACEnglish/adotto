Random annotation

tranno
```
chr1    248918871       248918909       17      2.3     80      1.6     CATATTTTCCTAATGAA
```

tie back to ucsd1
```
chr1    248918871       248918957       AAAATATGTTTATTAGG       17
```
These sequences don't match (though both are 17bp)

grch38 
```
>chr1:248918871-248918910
TCATATTTTCTTAATGAAAATATTTTCCTAATACACATAT
```

So I don't know what the ucsd1 sequence is
```
>chr1:248918871-248918957
TCATATTTTCTTAATGAAAATATTTTCCTAATACACATATATCAATGTGAGATTCATTTTTGTAAAAAAAATTATTTTTTTAATTTT

Is it the reverse compliment?
TCATATTTTCTTAATGAAAATATTTTCCTAATACACATAT
			  CCTAATAAACATATTTT <-- rc
		AAAATATGTTTATTAGG <-- orig
```
Nah, probably not

# Missing 1
region
```
chr10   27011197        27011258
```

ucsd1
```
chr10   27011222        27011233        AAAT    4
```

grch38
```
>chr10:27011197-27011258
ATGATACACAACAAAAACTAGTTGCTATTTATTTATTTTTGTTCAAACAAATGTTTTATGCT
			   AAAT <-- if its revcomp, maybe
```

pvcf
```
chr10   27011240        .       C       CA
```

