import re
import sys
import pysam

from math import ceil

def chunk_into_n(lst, n):
    size = ceil(len(lst) / n)
    return list(
        map(lambda x: lst[x * size:x * size + size],
            list(range(n)))
              )


f = pysam.FastaFile(sys.argv[1])
n_parts = 10
for pos, i in enumerate(chunk_into_n(f.references, n_parts)):
    with open(f'part{pos}.fasta', 'w') as fout:
        for ref in i:
            fout.write(f">{ref}\n")
            s = re.sub("(.{100})", "\\1\n", f[ref], 0, re.DOTALL).strip()
            fout.write(s + '\n')

