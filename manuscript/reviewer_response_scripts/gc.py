import pysam

f = pysam.FastaFile("/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa")
g = []
with open("pathoreg.txt", 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')
        chrom, start, end = data[:3]
        start = int(start)
        end = int(end)
        seq = f.fetch(chrom, start, end).upper()
        g.append((seq.count('G') + seq.count('C')) / len(seq))

import joblib
joblib.dump(g, 'gc.jl')
