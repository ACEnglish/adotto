#!/usr/bin/env python3
import sys
import json

suffix_name = '_' + sys.argv[1] if len(sys.argv) > 1 else ""
tot_cnt = 0
chr_rm_cnt = 0
sm_rm_cnt = 0
big_rm_cnt = 0
flip_rm_cnt = 0
span_rm = 0
span_kept = 0

known_chrs = [f"chr{n}" for n in range(1, 23)]
known_chrs.extend(["chrY", "chrX"])

for line in sys.stdin:
    tot_cnt += 1
    chrom, start, end  = line.strip().split('\t')
    start = int(start)
    end = int(end)
    span = abs(end - start)
    if chrom not in known_chrs:
        chr_rm_cnt += 1
        span_rm += span
        continue
    if end <= start:
        flip_rm_cnt += 1
        span_rm += span
        continue
    if abs(end - start) < 10:
        sm_rm_cnt += 1
        span_rm += span
        continue
    if span > 50000:
        big_rm_cnt += 1
        span_rm += span
        continue
    span_kept += span
    print(f"{chrom}\t{start}\t{end}")

span_total = span_kept + span_rm 
span_pct = span_rm / span_total * 100
rm_cnt = sm_rm_cnt + big_rm_cnt + flip_rm_cnt + chr_rm_cnt
pct_rm = rm_cnt / tot_cnt * 100

sys.stderr.write(f"""Count removed {rm_cnt} of {tot_cnt} = {tot_cnt - rm_cnt} ({pct_rm:.2f}%)
Span removed {span_rm} of {span_total} = {span_kept} ({span_pct:.2f}%)
Chromosome removed: {chr_rm_cnt}
Small removed: {sm_rm_cnt}
Big removed: {big_rm_cnt}
Coord remove: {flip_rm_cnt}
""")

out = {"total": tot_cnt,
       "removed": rm_cnt,
       "removed_pct": pct_rm,
       "chrom_rm": chr_rm_cnt,
       "small_rm": sm_rm_cnt,
       "big_rm": big_rm_cnt,
       "coord_rm": flip_rm_cnt,
       "span_total": span_total,
       "span_rm": span_rm,
       "span_kept": span_kept}
json.dump(out, open(f"merging_stats{suffix_name}.json", 'w'), indent=4)

