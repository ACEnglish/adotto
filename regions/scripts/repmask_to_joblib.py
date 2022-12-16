import sys

import joblib

from truvari.annotations.repmask import RepMask

from collections import defaultdict

lookup = defaultdict(list)
for i in sys.argv[1:]:

    data = RepMask.parse_output(i)
    for key in data:
        for item in data[key]:
            if item["RM_score"] >= 225:
                cls = item["RM_clsfam"].split('/')[0]
                lookup[key].append((item["RM_score"], cls))

joblib.dump(lookup, 'repmask_results.dict.jl')
