import sys
import pysam
import re


fasta = pysam.FastaFile(sys.argv[1])
chrom, rest = sys.argv[2].split(':')
start, end = rest.split('-')
start = int(start)
end = int(end)

oseq = fasta.fetch(chrom, start - 1, end)
oseq = re.sub("(.{60})", "\\1\n", oseq, 0, re.DOTALL)
print(f">ref_{chrom}_{start}_{end}\n{oseq}")
