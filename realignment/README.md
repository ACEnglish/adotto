In progress

bash msa_realign_pipe.sh reference.fa input.vcf.gz chrom:start-end

make a directory `msa_${chrom}:${start}-${end}`
Requires:
bcftools (with +mendelian)
vcftools (I'll package up vcf_compress)
samtools
mafft

python3
    pysam
    truvari

Procedure:
1) use bcftools consensus to create haplotypes
2) run mafft to perform multiple sequence alignment
    on unique haplotypes plus the reference sequence
3) turn msa result back into VCF

Limitations:
    Requires phased, sequence resolved variants
    Alters the variants such that counts won't be identical to inputs
        e.g. a caller has 10 variants, the truth set has 12, after msa_realign, there's 13 variants with 100% matching.
        will be hard to explain the increase in variants

    I think it might have more mendelian errors by variant count, but fewer errors by ALT bases
        (need to prove)

    errors will effect the variants
    It's hard to compare between runs. Imagine benchmaring callset1 to the benchmark and callset2 to the bencmark,
    because they may create different MSAs, different variant counts are produced. So if the two callsets have the same
    true-positive variants, and one FP each (that aren't the same), they could produce different counts of realigned
    variants. so precision/recall won't be comparable (e.g., if one of the FPs splits a TP into two but the other
    doesn't, then we got a precision1 = 1 / N and precision2 = 1 / (N + 1)
    
	Variants that overlap other variats aren't parsed by bcftools consensus.
	So gotta figure out how to warn about those
Positives:
    Alignming multiple samples at once helps tandem-repeats 'pile-up' better
        (need to prove)


Ideas:
    For the counting problem, what if the comparison is by bases.
    We could have true-negatives (ref homozygous). Presumably, there should be the same number of altenat

Remaining Items:
	The md5sum at the end needs to be identical. Otherwise I'm potentially altering the variants' sequences.
	Then I can pipeline all of this / wrap it up into a tool
	I want to log these stats
	Probably want a concatenate command somewhere
	And I need to figure out how to describe if the repeats are more unified
	Does this work for just a couple of haplotypes (as will be the case for benchmarking)
	Need to describe the process of
	1) bcftools merge (never norm, probably don't want multi-allelics, still
	2) run this script
	3) do the counting
	Shortcuts - only one variant, skip it
		  - only one (non-ref) haplotype, skip it
		  - minimum variant counts?
		  pre-anno of TR could help identify if there are any tandem-repeats and not just snps
	
	Maybe a number of entry_size * AF
	- That number should increase pretty good
	- and mean AF

	I kinda want to speed up msa2vcf.py - would be fun to make bindings of that
	Maybe a metric of AF std... like, I feel like the alleles should be come less 'jumpy' so shared-ness is more
	even, maybe?
	PCA? (very high level, probably already looks good.
	# count by shared.. hard

bedtools slop -b 500 -i <(grep -v _ ~/scratch/giabvntr/bed_files/baylor/merged.bed)  -g ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa.fai | bedtools merge -i - > simprep.bed
