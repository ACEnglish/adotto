import sys
import json
import pandas as pd

parts = []
for i in sys.argv[1:]:
    d = pd.Series(json.load(open(i, 'r')))
    d['name'] = i
    parts.append(d)

data = pd.concat(parts, axis=1).T
data.drop(columns=['gt_matrix'], inplace=True)
data.to_csv('/dev/stdout', sep='\t', header=True, index=False)
