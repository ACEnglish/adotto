import os
import pandas as pd
from io import StringIO
import glob

files = glob.glob("logs/cov.*.sh.log")
rows = []
for f in files:
    name = os.path.basename(f).split('.')[1:-2]
    str_data = StringIO()
    with open(f) as fh:
        for line in fh:
            if not line.startswith("#"):
                str_data.write(line)
    str_data.seek(0)
    m_dat = pd.read_csv(str_data, sep='\t')
    m_dat[["project", "sample", "sex"]] = name
    rows.append(m_dat)
data = pd.concat(rows)
data.to_csv("coverage_stats.txt", sep='\t', index=False)

