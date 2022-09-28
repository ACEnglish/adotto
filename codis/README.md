1) How many of the codis are already inside the tr catalog

 region counts:

 23 STR_Grch38_autosomes.bed
 29 STR_Grch38_sexchr.bed

 bedtools intersect -a STR_Grch38_autosomes.bed -b ~/scratch/code/adotto/regions/adotto_TRannotations_v0.2.bed.gz -c

 All regions are already inside the tr catalog

2) How many (and which) samples have adequate diploid coverage
Autosomes:
	sites with inadequate coverage
	| sites	|  samples |
	|-------|----------|
	|  0    |  78      |
	|  1    |  7       |
	|  2    |  1       |
The inadequately covered samples are:

	HG00732    2
	HG01258    1
	HG01361    1
	HG01978    1
	HG02622    1
	HG03486    1
	NA19240    1
	PGP1       1
SexChrs:
	sites with inadequate coverage
	| sites	|  samples |
	|-------|----------|
	| 0     | 8        |
	| 26    | 6        |
	| 24    | 6        |
	| 3     | 4        |
	| 25    | 3        |
	| 29    | 3        |
	| 1     | 3        |
	| 2     | 3        |
	| 22    | 1        |
	| 23    | 1        |
	| 27    | 1        |
	| 4     | 1        |
So only 8 samples have all adequately covered sites on the chrY codis sites

3) How many variants are inside the pVCF (per region)
	Only PentaE is without variants
	D1S1656	19
	D10S1248	3
	TH01	6
	vWA	14
	D12S391	24
	D13S317	5
	PentaE	0
	D16S539	3
	D18S51	12
	D19S433	13
	TPOX	2
	D2S441	9
	D2S1338	30
	D21S11	40
	PentaD	14
	D22S1045	6
	D3S1358	11
	FGA	21
	CSF1PO	7
	SE33	109
	D6S1043	7
	D7S820	13
	D8S1179	14

    chrY regions
    DYF387S1	7
    DYS19	36
    DYS385a	1
    DYS385b	1
    DYS389I	17
    DYS389II	47
    DYS390	6
    DYS391	2
    DYS392	2
    DYS393	2
    DYS437	6
    DYS438	4
    DYS439	69
    DYS448	13
    DYS456	1
    DYS458	9
    DYS460	0
    DYS461	5
    DYS481	2
    DYS505	27
    DYS522	0
    DYS533	1
    DYS549	3
    DYS570	6
    DYS576	5
    DYS612	9
    DYS635	8
    DYS643	2
    YGATAH4	2


Problems:

The coverage getter is only looking at the sites' start/end when it should be looking the whole thing.
This stuff needs to be cleaned up into a single report so we can present it more easily
