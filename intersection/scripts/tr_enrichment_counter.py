import sys
import truvari
import pysam
import pandas as pd
import numpy as np
import joblib

MISS_THRESHOLD = 0.10

regions = sys.argv[1] #pd.read_csv(sys.argv[1], sep='\t', header=None)
vcf = pysam.VariantFile(sys.argv[2])

m_regs = truvari.RegionVCFIterator(vcf, includebed=regions)

n_vals = np.zeros((4, 2))

for entry in vcf:
    #if entry.info["F_MISSING"] >= MISS_THRESHOLD:
    #    continue
    sz = truvari.entry_size(entry)
    row = 3
    if sz == 0:
        row = 0
    elif sz < 5:
        row = 1
    elif sz < 50:
        row = 2
    col = 0 if m_regs.include(entry) else 1
    n_vals[row, col] += 1
print(n_vals)
data = pd.DataFrame(n_vals, columns=["in", "out"], index=["SNP", "<5bp", "<50bp", "SV"])
joblib.dump(data, 'test.jl')

