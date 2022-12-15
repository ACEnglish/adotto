"""
Simplifies TR_region annotations and expands their columns

annotations are checked for 'nested', 'staggered', or 'isolated' categories
region boundaries are updated to try and ensure 25bp of buffer
new annotations are added to each region
- up_buff - how many bases upstream of the first annotation's start are known non-TR sequence
- dn_buff - how many bases downstream of the last annotation's end are known non-TR sequence
- n_filtered - how many annotations were removed from the region
- n_annos - how many annotations remain in the region
- n_subregions - how many subregions are in the region
- mu_purity - average purity of annotations in region
- pct_annotated - percent of the region's range (minus buffer) annotated
- ovl_flag - what annotation categories based on overlap filtering are inside the region


Additionally, each annotation has the following fields populated:
- purity - measured purity of the annotation over its span
- ovl_flag - a flag of the annotation's categories

the overlap flags are:
- iso [0x1] - if an annotation was isolated and by itself. True or False. May be indicative of a 'cleaner' repeat
- parent [0x2] - if annotation is part of a nested annotation filtering and is a parent repeat
- nested [0x4] - if annotation passed through nested annotation filtering and is a sub-repeat
- staggered [0x8] - if annotation passed through staggered annotation filtering. True or False. May be indicative of 'fuzzy'
  boundaries
"""
import sys
import json
import math
from functools import cmp_to_key

import pysam
import truvari
from intervaltree import IntervalTree

def iter_tr_regions(fn, region=None):
    """
    Read a repeats file with structure chrom, start, end, annotations.json
    returns generator of dicts
    """
    data = truvari.opt_gz_open(fn)
    for line in data:
        if isinstance(line, str):
            chrom, start, end, annos = line.strip().split('\t')
        else:
            chrom, start, end, annos = line
        start = int(start)
        end = int(end)
        annos = json.loads(annos)
        yield {"chrom": chrom,
               "start": start,
               "end": end,
               "annos": annos}


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
        anno['ovl_flag'] = 0
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
            anno['ovl_flag'] |= 0x4
            ret.append(anno)
    parent['ovl_flag'] |= 0x2
    ret.insert(0, parent)
    return ret

def viz_region(reg):
    """
    visualizes a region (poorly, but it's good enough)
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
    max_span['ovl_flag'] |= 0x8
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
        i['ovl_flag'] |= 0x1
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

BUFFER = 25
def get_bounds(annos):
    # Get the boundaries of the annotations
    start = min([_['start'] for _ in annos])
    end = max([_['end'] for _ in annos])
    return start, end

def annotate_purity(reg):
    seq = reference.fetch(reg['chrom'], reg['start'], reg['end'])
    scores = []
    for anno in reg["annos"]:
        s_start = anno['start'] - reg['start']
        s_end = anno['end'] - reg['start']
        sub_seq = seq[s_start: s_end]
        alt_seq = anno['repeat'] * int(math.floor(anno['copies']))
        remain = int((anno['copies'] % 1) * len(anno))
        alt_seq += anno['repeat'][:remain]
        if len(alt_seq) + len(sub_seq) == 0:
            print('what?')
            print(reg)
            print(anno)
            print(s_start, s_end)
            print(sub_seq)
            print(alt_seq1)
            print(alt_seq)
            
            assert False
        score = truvari.seqsim(sub_seq, alt_seq) * 100
        anno['purity'] = round(score)
        scores.append(score)
    return round(sum(scores) / len(scores))

def pct_annotated(reg, start, end):
    m_tree = IntervalTree()
    for i in reg['annos']:
        m_tree.addi(i['start'], i['end'])
    m_tree.merge_overlaps()
    tot_anno = sum([_.end - _.begin for _ in m_tree])
    return round(tot_anno / (end - start) * 100)

def write_region(reg):
    out_str = []
    for key in ['chrom', 'start', 'end', 'ovl_flag', 'up_buff', 'end_buff', 'n_filtered', \
                'n_annos', 'n_subregions', 'mu_purity', 'pct_annotated']:
        out_str.append(str(reg[key]))
    out_str.append(json.dumps(reg['annos']))
    print("\t".join(out_str))
        
def get_blur(reg, new_start, new_end):
    """
    Returns the percent of buffer bases known to not contain TR sequence
    """
    tree = IntervalTree()
    for anno in reg["annos"]:
        tree.addi(anno['start'], anno['end'])
    tree.merge_overlaps()
    up_cnt = 0
    for i in range(max(reg['start'], new_start), new_start + BUFFER):
        up_cnt += int(tree.overlaps(i))
    dn_cnt = 0
    for i in range(new_end - BUFFER, min(new_end, reg['end'])):
        dn_cnt += int(tree.overlaps(i))
    up_cnt = round((1 - (up_cnt / BUFFER)) * 100)
    dn_cnt = round((1 - (dn_cnt / BUFFER)) * 100)
    return up_cnt, dn_cnt

def get_buff(reg, f_start, l_end):
    """
    Try to put at least 25bp of buffer on the ends.
    """
    tree = IntervalTree()
    for anno in reg["annos"]:
        tree.addi(anno['start'], anno['end'])
    tree.merge_overlaps()
    start = f_start
    up_buff = 0
    best_upstream = max(f_start - BUFFER, reg['start']) # the furthest we've checked
    while not tree.overlaps(start - 1) and start > best_upstream:
        start -= 1
        up_buff += 1
    end = l_end
    dn_buff = 0
    best_dnstream = min(l_end + BUFFER, reg['end'])
    while not tree.overlaps(end) and end < best_dnstream:
        end += 1
        dn_buff += 1
    return start, end, up_buff, dn_buff

if __name__ == '__main__':
    in_anno, in_ref = sys.argv[1:]
    reference = pysam.FastaFile(in_ref)
    removed = 0
    total = 0
    for reg in iter_tr_regions("adotto_TRannotations_v0.3.bed.gz"):
        total += 1
        # Filter homopolymers
        in_anno_cnt = len(reg['annos'])
        out_annos = simplify_region([_ for _ in reg['annos'] if len(_['repeat']) != 1])
        if len(out_annos) == 0:
            removed += 1
            continue
        reg['ovl_flag'] = 0
        for i in out_annos:
            reg['ovl_flag'] |= i['ovl_flag']
        #print(reg)
        orig_start, orig_end = reg['start'], reg['end']
        first_start, last_end = get_bounds(out_annos)
        new_start, new_end, s_buff, e_buff = get_buff(reg, first_start, last_end)
        reg['up_buff'], reg['dn_buff'] = s_buff, e_buff
        reg['start'], reg['end'] = new_start, new_end
        reg['annos'] = out_annos
        reg['n_filtered'] = in_anno_cnt - len(out_annos)
        reg['n_annos'] = len(out_annos)
        reg['n_subregions'] = len([_ for _ in get_subregions(out_annos)])
        reg['mu_purity'] = annotate_purity(reg)
        reg['pct_annotated'] = pct_annotated(reg, first_start, last_end)
        
        write_region(reg)
    sys.stderr.write(f"removed {removed} regions from {total}\n")
