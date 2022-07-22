#!/usr/bin/env python3
import sys
sm_rm_cnt = 0
big_rm_cnt = 0
flip_rm_cnt = 0
tot_cnt = 0
for line in sys.stdin:
    tot_cnt += 1
    chrom, start, end  = line.strip().split('\t')
    start = int(start)
    end = int(end)
    span = abs(end - start)
    if abs(end - start) < 10:
        sm_rm_cnt += 1
        continue
    if span > 50000:
        big_rm_cnt += 1
        continue
    if end <= start:
        flip_rm_cnt += 1
        continue

    print(f"{chrom}\t{start}\t{end}")

rm_cnt = sm_rm_cnt + big_rm_cnt + flip_rm_cnt
pct_rm = rm_cnt / tot_cnt * 100
out_log = f"""Removed {rm_cnt} of {tot_cnt} = {tot_cnt - rm_cnt} ({pct_rm:.2f}%)
Small removed: {sm_rm_cnt}
Big removed: {big_rm_cnt}
Coord remove: {flip_rm_cnt}
"""

with open("merging_stats.txt", 'w') as fout:
    fout.write(out_log)
sys.stderr.write(out_log)

