"""
Consolidates the stats generated by the merges into a reasonable table.txt
"""
import os
import re
import sys
import json
import pandas as pd
in_dir = sys.argv[1] if len(sys.argv) > 1 else 'data/'

def spans_sum_read(fn):
    ret = {}
    with open(fn) as fh:
        for line in fh:
            name, val = re.split('\s+', line.strip())
            ret[name] = int(val)
    return pd.Series(ret)

def merge_json_read(fn):
    return pd.Series(json.load(open(fn)))
"""
data/
    input_spans.txt
    input_spans_merged.txt
    input_spans_slop.txt
    merging_stats.json
    merging_stats_slop.json
    [baylor|giab|pacbio|ucsd1|ucsd2]
        input_spans.txt
        merging_stats.json
"""

rows = []
for source in ['baylor', 'giab', 'pacbio', 'ucsd1', 'ucsd2']:
    spans = spans_sum_read(os.path.join(in_dir, source, 'input_spans.txt'))
    merge = merge_json_read(os.path.join(in_dir, source, 'merging_stats.json'))
    my_row = pd.DataFrame(pd.concat([spans, merge])).T
    my_row["source"] = source
    rows.append(my_row)

spans = spans_sum_read(os.path.join(in_dir, "input_spans.txt"))
merge = merge_json_read(os.path.join(in_dir, 'merging_stats.json'))
my_row = pd.DataFrame(pd.concat([spans, merge])).T
my_row["source"] = 'grand'
rows.append(my_row)

spans = spans_sum_read(os.path.join(in_dir, "input_spans_merged.txt"))
merge = merge_json_read(os.path.join(in_dir, 'merging_stats_slop.json'))
my_row = pd.DataFrame(pd.concat([spans, merge])).T
my_row["source"] = 'grand slop'
rows.append(my_row)


spans = spans_sum_read(os.path.join(in_dir, "input_spans_slop.txt"))
my_row = pd.DataFrame([spans])
my_row["source"] = "pre-gap"
rows.append(my_row)

spans = spans_sum_read(os.path.join(in_dir, "input_spans_final.txt"))
my_row = pd.DataFrame([spans])
my_row["source"] = "final"
rows.append(my_row)
data = pd.concat(rows)
data.to_csv("/dev/stdout", sep='\t', index=False)

