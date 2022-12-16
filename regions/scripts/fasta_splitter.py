import re
import pysam

f = pysam.FastaFile("tr_regions.v0.3.fasta")


first = f.references[:453472]
second = f.references[453472:453472 + 453472]
third = f.references[-453472:]

with open('first.fasta', 'w') as fout:
    for ref in first:
        s = re.sub("(.{100})", "\\1\n", f[ref], 0, re.DOTALL)
        fout.write(f">{ref}\n{s}\n")

with open('second.fasta', 'w') as fout:
    for ref in second:
        s = re.sub("(.{100})", "\\1\n", f[ref], 0, re.DOTALL)
        fout.write(f">{ref}\n{s}\n")

with open('third.fasta', 'w') as fout:
    for ref in third:
        s = re.sub("(.{100})", "\\1\n", f[ref], 0, re.DOTALL)
        fout.write(f">{ref}\n{s}\n")
