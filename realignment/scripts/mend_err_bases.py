import sys
import pysam
import truvari


v = pysam.VariantFile(sys.stdin)
tot_bases = 0
err_bases = 0
err_vars = 0
entry_count = 0
for entry in v:
    sz = truvari.entry_size(entry)
    entry_count += 1
    if entry.info["MERR"] != 0:
        err_bases += truvari.entry_size(entry)
        err_vars += 1
    tot_bases += sz
print("total variants:", entry_count)
print("total bases:", tot_bases)
print("menderr bases:", err_bases)
print("menderr vars:", err_vars)
        
