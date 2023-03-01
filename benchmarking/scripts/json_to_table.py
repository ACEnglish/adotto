import os
import sys
import json
import pandas as pd

rows = []
for i in sys.argv[1:]:
    name, fn = i.split(':')
    m_data = pd.DataFrame([json.load(open(fn))])
    m_data['program'] = name
    m_data['method'] = 'refine' if 'refine' in os.path.basename(fn) else 'bench'
    rows.append(m_data)

data = pd.concat(rows, axis=0, ignore_index=True)
data.drop(columns=["gt_matrix"], inplace=True)
data.to_csv("/dev/stdout", sep='\t', index=False)
