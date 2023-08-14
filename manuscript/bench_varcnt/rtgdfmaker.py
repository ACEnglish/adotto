import os
import sys
import joblib
import truvari
import pandas as pd

dir_name = sys.argv[1]

files = [('fn', 'fn.vcf.gz'), 
         ('fp', 'fp.vcf.gz'),
         ('tp-base', 'tp-baseline.vcf.gz'),
         ('tp', 'tp.vcf.gz')]

parts =  []
for state, fn in files:
    data = truvari.vcf_to_df(os.path.join(dir_name, fn), True, True)
    data['state'] = state
    parts.append(data)

data = pd.concat(parts)
joblib.dump(data, os.path.join(dir_name, 'data.jl'))
