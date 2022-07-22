import sys
import truvari
import pandas as pd

from truvari.vcf2df import pull_samples

d = truvari.vcf_to_df(sys.argv[1], with_fmt=True, sample=pull_samples(sys.argv[1]))

gt_cols = [_ for _ in d.columns if _.endswith("_GT")]

rows = []
for i in gt_cols:
    d[i] = d[i].apply(lambda x: truvari.get_gt(x).name)
    x = d[i].value_counts()
    he = x.loc['HET'] if 'HET' in x else 0
    ho = x.loc['HOM'] if 'HOM' in x else 0
    x['tot'] = he + ho
    x['samp'] = i
    rows.append(x)

print(pd.concat(rows, axis=1).loc['tot'].astype(int).describe())

