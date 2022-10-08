import os
import sys
import glob
import truvari
import joblib
import pandas as pd

in_folder = sys.argv[1]
out_file = sys.argv[2]

outputs = glob.glob(os.path.join(in_folder, "*", "output.vcf.gz"))
#outputs = ["results/cmrg_phab/chr10:101354033-101366711/output.vcf.gz"]
parts = []
summary = []
for i in outputs:

    data = truvari.vcf_to_df(i, with_fmt=True, sample=["HG002", "pHG002"])
    bpres = data['HG002_GT'].apply(lambda x: truvari.get_gt(x) in [truvari.GT.HET, truvari.GT.HOM])
    cpres = data['pHG002_GT'].apply(lambda x: truvari.get_gt(x) in [truvari.GT.HET, truvari.GT.HOM])
    data['state'] = None
    data.loc[bpres == cpres, 'state'] = 'tp'
    data.loc[(bpres != cpres) & bpres, 'state'] = 'fn'
    data.loc[(bpres != cpres) & cpres, 'state'] = 'fp'

    cnts = data.groupby(['state']).size()
    region = i.split('/')[-2]
    tp = cnts['tp'] if 'tp' in cnts else 0
    fp = cnts['fp'] if 'fp' in cnts else 0
    fn = cnts['fn'] if 'fn' in cnts else 0
    summary.append([region, tp, fp, fn])
    parts.append(data)

data = pd.concat(parts)
joblib.dump(data, out_file)


d = pd.DataFrame(summary, columns=['region', 'tp', 'fp', 'fn'])
d['precision'] = d['tp'] / (d['tp'] + d['fp'])
d['recall'] = d['tp'] / (d['tp'] + d['fn'])
d['f1'] = 2 * (d['recall'] * d['precision']) / (d['recall'] + d['precision'])
joblib.dump(d, out_file + '.summary.jl')
