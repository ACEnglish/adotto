"""
Simplifies TR_region annotations

annotations are checked for 'nested', 'staggered', or 'isolated' conditions
region boundaries are updated to ensure 25bp of buffer
new annotations are added to each region
- n_filtered - how many annotations were removed from the region
- n_annos - how many annotations remain in the region
- n_subregions - how many subregions are in the region
- mu_purity - average purity of annotations in region
- pct_annotated - percent of the region's range (minus buffer) annotated

Additionally, each annotation has the following fields populated:
- purity - measured purity of the annotation over its span
- nested - if annotation passed through nested annotation filtering. 0 is a parent repeat. 1 is a sub-repeat
- staggered - if annotation passed through staggered annotation filtering. True or False. May be indicative of 'fuzzy'
  boundaries
- iso - if an annotation was isolated and by itself. True or False. May be indicative of a 'cleaner' repeat
"""
import json
import math
import sys
from functools import cmp_to_key

from truvari.annotations.trf import iter_tr_regions
from suffix_tree import Tree
from intervaltree import IntervalTree
import numpy as np

"""
1. for each region
2. :: remove homopolymer only annotations ::

2. Separate out sub-regions (collect stats)
3. for each sub-region
    
    are the annotations 'pure-tree(?)' - where every interval is a 
    sub-interval of some biggest span?
    if yes:
        just do common_substring filtering - if the motif isn't part of the common_substring,
        remove it
    if no:
        I haven't seen this case, yet
"""
def same_intv(intv, tree):
    """
    Return if intv is in tree (no data comparison
    """
    for i in tree:
        if i.begin == intv.begin and i.end == intv.end:
            return True
    return False

def get_subregions(annos):
    """
    For a set of ranges, separate them into distinct, non-overlapping regions
    # also need a way to pull out how many we're working with?
    """
    i_tree = IntervalTree()
    for pos, anno in enumerate(annos):
        i_tree.addi(anno['start'], anno['end'], pos)
    
    m_tree = i_tree.copy()
    m_tree.merge_overlaps()
    for i in m_tree:
        m_y = [annos[_.data] for _ in i_tree.overlap(i)]
        # Single annotation
        if len(m_y) == 1:
            m_t = 'isolated'
        elif same_intv(i, i_tree): # merged interval is an exact interval of annos
            m_t = 'nested'
        else:
            m_t = 'staggered'
        yield m_y, m_t

def build_common_tree(all_strs):
    """
    Builds dict with key: repeat value: is part of the common substrings of the suffix tree
    """
    m_tree = Tree()

    for i, j in enumerate(all_strs):
        m_tree.add(i, j * 2)
    common = []
    #for C, path in sorted(m_tree.maximal_repeats()):
    for k, length, path in m_tree.common_substrings():
        cur_row = []
        x = "".join([str(_) for _ in path])
        for i in all_strs:
            cur_row.append(i in x)
        common.append(cur_row)
    common = np.array(common)
    is_keep = dict(zip(all_strs, common.any(axis=0)))
    return is_keep


def anno_sorter(a1, a2):
    if a1[0] > a2[0]:
        return 1
    if a1[0] < a2[0]:
        return -1
    # if the boundaries are identical, I want the smaller motif?
    #if a1[1]['start'] == a2[1]['start'] and a1[1]['end'] == a2[1]['end']:
    #    return len(a2[1]['repeat']) 
    return a1[1]["score"] - a2[1]["score"]

def make_span_sorted_list(annos):
    """
    Sort the annotations by their span then score, return list of tuples
    """
    srt = []
    for anno in annos:
        srt.append((anno['end'] - anno['start'], anno))
    # in the case of a span tie, take the one with the higher score
    srt.sort(key=cmp_to_key(anno_sorter))
    return srt

def simplify_nested(annos):
    """
    Get the 'parent' annotation (assume longest for now)
    Add that to the tree
    is_keep will be tree.find(i)
    """
    srt = make_span_sorted_list(annos)
    
    parent = srt[-1][1] # assuming longest span is the parent.
    m_tree = IntervalTree()
    pos = parent['start']
    for i in range(int(parent['copies'])):
        end = pos + int(parent['period'])
        m_tree.addi(pos, end)
        pos = end
    ret = []
    for l, anno in srt[:-1]:
        if len(m_tree.overlap(anno['start'], anno['end'])) == 1:
            anno['nested'] = 1
            ret.append(anno)
    parent['nested'] = 0
    ret.insert(0, parent)
    return ret

def viz_region(reg):
    """
    visualizes a region
    """
    #print(">chr1:52111-52334")
    #print("ACATGATTTTTTTCTTTGCTGTTCTTGTCTAATTGTTATTAATAATTAATAAATAACTTATGATCTAATTGTTATTAATAATAACTTATCATCACATGATTTATTAATAAATTAATAAATAACTTATTATCACCGCATTTCCCCAATTCATTTATCTTTCTTTCATTTTCTCTCTTTGTGTGTTTTCTGTCTTCATATTTCAGCACTTGCCACATATTTCCCAC")
    start_pos = reg['start']
    for pos, anno in enumerate(reg['annos']):
        buffer = " " * (anno['start'] - start_pos)
        # I want to alternate upper, lower
        seq = ""
        mode = True
        for i in range(math.floor(anno['copies'])):
            if mode:
                seq += anno['repeat']
            else:
                #seq += "[cyan]" + anno['repeat'].lower() + "[/]"
                seq += anno['repeat'].lower()
            mode = not mode
    
        anno['end'] - anno['start']
        remain = anno['end'] - anno['start']
        remain %= int(anno['period'])
        #if remain != 0:
        #seq = seq[:-remain]
        print(buffer + seq, '(', anno['start'] - start_pos, anno['end'] - start_pos, ')', anno['score'])

def simplify_staggered(annos):
    """
    Return a list of them..
    """
    srt = make_span_sorted_list(annos)
    max_span = srt[-1][1]
    max_span['stag'] = True
    keeps = [max_span]
    for l, intv in srt[:-1]:
        if max_span["start"] <= intv["start"] and intv["end"] <= max_span["end"]:
            keeps.append(intv)
        elif intv["start"] > max_span["end"] or intv["end"] < max_span["start"]:
            keeps.append(intv)
        # filtered
    for i in simplify_region(keeps):
        yield i

def annotate_iso(annos):
    for i in annos:
        i['iso'] = True
    return annos

def simplify_region(annos):
    """
    Parent procedure
    Returns list of new annotations
    """
    ret = []
    for sub_region_annos, tr_type in get_subregions(annos):
        if tr_type == "isolated":
            ret.extend(annotate_iso(sub_region_annos))
        elif tr_type == 'nested':
            ret.extend(simplify_nested(sub_region_annos))
        elif tr_type == 'staggered':
            ret.extend(simplify_staggered(sub_region_annos))
    return ret

def get_bounds(annos):
    # Get the boundaries of the annotations
    BUFFER = 25
    start = min([_['start'] for _ in annos]) - BUFFER
    end = max([_['end'] for _ in annos]) + BUFFER
    return start, end

def annotate_purity(annos):
    # TODO - I have this code somewhere
    return 1

def pct_annotated(reg):
    m_tree = IntervalTree()
    for i in reg['annos']:
        m_tree.addi(i['start'], i['end'])
    m_tree.merge_overlaps()
    tot_anno = sum([_.end - _.begin for _ in m_tree])
    return round(tot_anno / (reg['end'] - reg['start'] - 25), 3)

def write_region(reg):
    out_str = []
    for key in ['chrom', 'start', 'end', 'n_filtered', 'n_annos', 'n_subregions', 'mu_purity', 'pct_annotated']:
        out_str.append(str(reg[key]))
    out_str.append(json.dumps(reg['annos']))
    print("\t".join(out_str))
        

if __name__ == '__main__':
    # need to do 
    removed = 0
    total = 0
    for reg in iter_tr_regions("adotto_TRannotations_v0.3.bed.gz"):
        total += 1
        #viz_region(reg)

        # Filter homopolymers
        in_anno_cnt = len(reg['annos'])
        out_annos = simplify_region([_ for _ in reg['annos'] if len(_['repeat']) != 1])
        if len(out_annos) == 0:
            removed += 1
            continue
        reg['start'], reg['end'] = get_bounds(out_annos)
        reg['n_filtered'] = in_anno_cnt - len(out_annos)
        reg['n_annos'] = len(out_annos)
        reg['n_subregions'] = len([_ for _ in get_subregions(out_annos)])
        reg['mu_purity'] = annotate_purity(out_annos)
        reg['annos'] = out_annos
        reg['pct_annotated'] = pct_annotated(reg)
        
        write_region(reg)
    sys.stderr.write(f"removed {removed} regions from {total}\n")
