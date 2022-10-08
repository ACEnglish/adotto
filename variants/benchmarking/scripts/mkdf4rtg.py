import joblib
import pandas as pd
import sys
import os
import truvari

in_folder = sys.argv[1]
files = [("fn.vcf.gz", "fn"), ("fp.vcf.gz", "fp"), ("tp-baseline.vcf.gz", "tpbase"), ("tp.vcf.gz", "tp")]

parts = []
for fn, state in files:
    fn = os.path.join(in_folder, fn)
    d = truvari.vcf_to_df(fn, True)
    d['state'] = state
    parts.append(d)

data = pd.concat(parts)
joblib.dump(data, os.path.join(in_folder, "data.jl"))
