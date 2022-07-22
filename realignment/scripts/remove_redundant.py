import sys
import pysam
from collections import defaultdict
f = pysam.FastaFile(sys.argv[1])
seqs = defaultdict(list)
total = 0
redundant = 0
kept = 0
for i in f.references:
    seq = f.fetch(i)
    total += 1
    if seq in seqs:
        redundant += 1
    else:
        kept += 1
    seqs[seq].append(i)

for seq in seqs:
    key = ";".join(seqs[seq])
    print(f">{key}\n{seq}")

sys.stderr.write(f"{redundant} redundant. {kept} kept. {total} total\n")
