import sys
import json
import pysam
import truvari
from collections import Counter
from truvari.annotations.lcr import sequence_entropy
def filter_TR_entry(entry):
    """
    Returns True if this should be filtered
    """
    if "TRFdiff" not in entry.info or entry.info["TRFdiff"] == None:
        return True
    if entry.info["TRFperiod"] == 1:
        return True
    if len(entry.alts[0]) > len(entry.ref):
        seq = entry.alts[0][1:]
    else:
        seq = entry.ref[1:]
    if sequence_entropy(seq) < 0.25:
        return True
    return False

if __name__ == '__main__':
    v = pysam.VariantFile(sys.argv[1])
    #bed = sys.argv[2]
    header = v.header.copy()
    header.add_line('##FILTER=<ID=NONTR,Description="NonTR variant">')
    o = pysam.VariantFile('/dev/stdout', 'w', header=header)
    for entry in v:
        if truvari.entry_size(entry) < 5:
            o.write(entry)
            continue
        entry.translate(header)
        if filter_TR_entry(entry):
            entry.filter.clear()
            entry.filter.add('NONTR')
        o.write(entry)
