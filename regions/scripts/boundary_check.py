"""
Checks the boundaries of regions
expect 25bp of non-TR sequence buffering all annotations
"""
import sys
from truvari.annotations.trf import iter_tr_regions
from collections import Counter
in_bed = sys.argv[1]
cnt = Counter()
thresh = 25
for region in iter_tr_regions(in_bed):
    has_prob = False
    for anno in region["annos"]:
        if abs(anno["start"] - region["start"]) < thresh or abs(anno["end"] - region["end"]) < thresh:
            has_prob = True
    cnt[has_prob] += 1
import json

print(json.dumps(cnt, indent=4))


