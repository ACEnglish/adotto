import pysam
import sys

f = pysam.FastaFile(sys.argv[1])

for i in f.references:
    print(f">{i}\n{f[i].upper()}")
