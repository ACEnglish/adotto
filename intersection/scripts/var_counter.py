import sys
import truvari
import pysam
import pandas as pd

MISS_THRESHOLD = 0.10
regions = pd.read_csv(sys.argv[1], sep='\t', header=None)
vcf = pysam.VariantFile(sys.argv[2])

n_rows = []
#cnt = 10
for idx, region in regions.iterrows():
    lt5 = 0
    ge5 = 0
    for entry in vcf.fetch(region[0], region[1], region[2]):
        if entry.info["F_MISSING"] >= MISS_THRESHOLD:
            continue
        ent_start, ent_end = truvari.entry_boundaries(entry)
        if region[1] <= ent_start and ent_end <= region[2]:
            sz = truvari.entry_size(entry)
            if sz < 5:
                lt5 += 1
            else:
                ge5 += 1
    n_rows.append([lt5, ge5])
    #if len(n_rows) >= cnt:
    #    break
n_rows = pd.DataFrame(n_rows, columns=["var_lt5", "var_ge5"], index=regions.index)
regions = regions.join(n_rows)
regions.to_csv("/dev/stdout", sep='\t', header=False, index=False)


