import os
import json
import pandas as pd
import joblib

bench_parts = {}
refine_var = {}
refine_reg = {}
with open("paths.txt", 'r') as fh:
    fh.readline()
    for line in fh:
        data = line.strip().split('\t')
        p = os.path.join(data[1], 'summary.json')
        bench_parts[data[0]] = json.load(open(p))
        p = os.path.join(data[1], 'refine.variant_summary.json')
        refine_var[data[0]] = json.load(open(p))
        p = os.path.join(data[1], 'refine.region_summary.json')
        refine_reg[data[0]] = json.load(open(p))

data = {'bench': pd.DataFrame(bench_parts).T,
        'refine_v': pd.DataFrame(refine_var).T,
        'refine_r': pd.DataFrame(refine_reg).T
        }
joblib.dump(data, 'all.summary.jl')
