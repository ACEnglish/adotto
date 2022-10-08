import truvari
import pandas as pd
import joblib


data = truvari.vcf_to_df('/dev/stdin', with_fmt=True, sample=["HG002", "pHG002"])
bpres = data['HG002_GT'].apply(lambda x: truvari.get_gt(x) in [truvari.GT.HET, truvari.GT.HOM])
cpres = data['pHG002_GT'].apply(lambda x: truvari.get_gt(x) in [truvari.GT.HET, truvari.GT.HOM])
data['state'] = None
data.loc[bpres == cpres, 'state'] = 'tp'
data.loc[(bpres != cpres) & bpres, 'state'] = 'fn'
data.loc[(bpres != cpres) & cpres, 'state'] = 'fp'

cnts = data.groupby(['state']).size()
tp = cnts['tp'] if 'tp' in cnts else 0
fp = cnts['fp'] if 'fp' in cnts else 0
fn = cnts['fn'] if 'fn' in cnts else 0

recall = tp / (tp + fn)
precision = tp / (tp + fp)
f1 = 2 * (recall * precision) / (recall + precision)
summary = pd.Series([tp, fp, fn, recall, precision, f1], index=['tp', 'fp', 'fn', 'recall', 'precision', 'f1'])

joblib.dump(data, 'phab.raw.jl')
print(summary)



