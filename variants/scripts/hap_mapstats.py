"""
Given an adotto output directory with alignment/variant calling results from `map_haplo.sh`
Consolidate stats into a pandas Series
"""
import os
import sys
import joblib
import pandas as pd
parts = []
for in_dir in sys.argv[1:]:
    data = {}
    name = os.path.basename(in_dir.rstrip('/'))
    with open(os.path.join(in_dir, 'aln.paf.stats.txt'), 'r') as fh:
        for line in fh:
            key, val = line.strip().split(':')
            data[key.strip()] = int(val.strip())

    data = pd.Series(data)
    cov = pd.read_csv(os.path.join(in_dir, "cov.bed"), sep='\t', names=["chrom", "start", "end", "coverage"])
    cov['span'] = cov['end'] - cov['start']
    summary_coverage = cov.groupby('coverage')['span'].sum()

    data['reference bases single covered'] = summary_coverage.loc[1]
    data['reference bases uncovered']  = summary_coverage.loc[0]
    data['reference bases multi covered'] = summary_coverage[~summary_coverage.index.isin([0, 1])].sum() 
    data.name = name
    parts.append(data)

data = pd.DataFrame(parts)
joblib.dump(data, 'hap_mapstats.jl')
