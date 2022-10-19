import sys
from truvari.annotations.trf import iter_tr_regions
from intervaltree import IntervalTree

def tree_and_merge(line):
    
    m_tree = IntervalTree()
    for i in line['annos']:
        m_tree.addi(i['start'], i['end'])
    line['regcnt'] = len(m_tree)
    m_tree.merge_overlaps()
    line['mrgcnt'] = len(m_tree)
    return line

parts = []
for entry in iter_tr_regions(sys.argv[1]):
    entry = tree_and_merge(entry)
    parts.append(entry)
    #entry['mrgcnt'])
import joblib
joblib.dump(parts, 'tr_regmrgcnts.jl')
