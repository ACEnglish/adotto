"""
Simplifies TR_region annotations and expands their columns

annotations are checked for 'nested', 'staggered', or 'isolated' categories
region boundaries are updated to try and ensure 25bp of buffer
new annotations are added to each region

- ovl_flag - what annotation categories based on overlap filtering are inside the region
- up_buff - how many bases upstream of the first annotation's start are known non-TR sequence
- dn_buff - how many bases downstream of the last annotation's end are known non-TR sequence
- hom_span - how many bases of the region were found to be homopolymer repeats
- n_filtered - how many annotations were removed from the region
- n_annos - how many annotations remain in the region
- n_subregions - how many subregions are in the region
- mu_purity - average purity of annotations in region
- pct_annotated - percent of the region's range (minus buffer) annotated
- intersersed - name of interspersed repeat found within region by RepeatMasker
- patho - name of gene affected by a pathogenic tandem repeat within region
- codis - name of codis site contained within region
(later adds
    gene_flag
    biotype
)
then there's annotations json
the overlap flags are:
- iso [1] - if an annotation was isolated and by itself. True or False. May be indicative of a 'cleaner' repeat
- parent [2] - if annotation is part of a nested annotation filtering and is a parent repeat
- nested [4] - if annotation passed through nested annotation filtering and is a sub-repeat
- staggered_dn [8] - if annotation passed through staggered annotation filtering on its downstream boundary
- staggered_up [16] - if annotation passed through staggered annotation filtering on its upstream boundary
"""
import sys
import json
import math
from functools import cmp_to_key

import pysam
import joblib
import truvari
import pandas as pd
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
HEADER = ['chrom', 'start', 'end', 'ovl_flag', 'up_buff', 'dn_buff', 'hom_span', 'n_filtered', 'n_annos', 'n_subregions',
          'mu_purity', 'pct_annotated', 'interspersed', 'patho', 'codis', 'annos']
def iter_tr_regions_1(fn):
    """
    Reads an already processed repeats file
    """
    data = truvari.opt_gz_open(fn)
    for line in data:
        line = {x:y for x,y in zip(HEADER, line.strip().split('\t'))}
        line['annos'] = json.loads(line['annos'])
        for i in HEADER[1:-4]:
            line[i] = int(line[i])
        yield line

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
        if 'ovl_flag' not in anno:
            anno['ovl_flag'] = 0
        if 'repeat' in anno:
            anno['motif'] = anno['repeat']
            del anno['repeat']
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
    # if the boundaries are identical, I want the higher score
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
            anno['ovl_flag'] |= 2 # nested
            ret.append(anno)
    parent['ovl_flag'] |= 4 # parent
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
                seq += anno['repeat'] if 'repeat' in anno else anno['motif']
            else:
                #seq += "[cyan]" + anno['repeat'].lower() + "[/]"
                seq += anno['repeat'].lower() if 'repeat' in anno else anno['motif'].lower()
            mode = not mode
        
        if mode:
            seq += seq[:(max(0, anno['end'] - anno['start'] - len(seq)))]
        else:
            seq += seq[:(max(0, anno['end'] - anno['start'] - len(seq)))].lower()

        #anno['end'] - anno['start']
        #remain = anno['end'] - anno['start']
        #remain %= int(anno['period'])
        #if remain != 0:
            #seq = seq[:-remain]
        print(buffer + seq, '(', anno['start'] - start_pos, anno['end'] - start_pos, ')', anno['score'])

def overlap_amount(astart, aend, bstart, bend):
    """
    Calculates the number of bases shared between two ranges
    :param `astart`: First range's start position
    :type `astart`: int
    :param `aend`: First range's end position
    :type `aend`: int
    :param `bstart`: Second range's start position
    :type `bstart`: int
    :param `bend`: Second range's end position
    :type `bend`: int

    :return: overlap amount
    :rtype: int
    """
    if astart >= bstart and aend <= bend:
        return 1
    ovl_start = max(astart, bstart)
    ovl_end = min(aend, bend)
    if ovl_start < ovl_end:  # Otherwise, they're not overlapping
        ovl_amt = ovl_end - ovl_start
    else:
        ovl_amt = 0
    return ovl_amt

def staggered_chooser(base, annos, up=True):
    """
    From the set of annos which would extend out the furthest. 
    Ties broken by having the least overlap with the base_anno
    """
    def up_sorter(anno_a, anno_b):
        # Pos - a goes after b
        # Neg - b goes after a
        if anno_a['start'] == anno_b['start']:
            a_ovl = overlap_amount(base['start'], base['end'], anno_a['start'], anno_a['end'])
            b_ovl = overlap_amount(base['start'], base['end'], anno_b['start'], anno_b['end'])
            if a_ovl == b_ovl:
                # score is final decider
                return anno_b['score'] - anno_a['score']
            # overlap is second decider
            return b_ovl - a_ovl
        # furthest up stream is first decider
        return anno_a['start'] - anno_b['start']

    def dn_sorter(anno_a, anno_b):
        # Pos - a goes after b
        # Neg - b goes after a
        if anno_a['end'] == anno_b['end']:
            a_ovl = overlap_amount(base['start'], base['end'], anno_a['start'], anno_a['end'])
            b_ovl = overlap_amount(base['start'], base['end'], anno_b['start'], anno_b['end'])
            if a_ovl == b_ovl:
                # score is final decider
                return anno_b['score'] - anno_a['score']
            # overlap is second decider
            return b_ovl - a_ovl
        # furthest dn stream is first decider
        return - (anno_a['end'] - anno_b['end'])

    m_srt = up_sorter if up else dn_sorter
    annos.sort(key=cmp_to_key(m_srt))
    # the first is the one we want
    return annos[0]

def mutual_contain(anno_a, anno_b):
    """
    returns True if anno_a (a tuple of coords) contains anno_b or vice versa
    """
    if anno_a[0] < anno_b[0] < anno_b[1] < anno_a[1]:
        return True
    if anno_b[0] < anno_a[0] < anno_a[1] < anno_b[1]:
        return True
    return False

def simplify_staggered(annos):
    """
    Return a list of them..
    # I need to separate them so that they exhaust.
    So upstream, downstream, and nested
    # all anchored around the max_span
    """
    srt = make_span_sorted_list(annos)
    max_span = srt[-1][1]
    keep_center = [max_span]
    up_mspan_copy = max_span["start"], max_span["start"] + max_span["period"]
    dn_mspan_copy = max_span["end"] - max_span["period"], max_span["end"]
    keep_up = []
    keep_dn = []
    up_need_to_choose = [] 
    dn_need_to_choose = [] 
    for l, intv in srt[:-1]:
        intv["ovl_flag"] |= 8
        # contained
        if max_span["start"] <= intv["start"] and intv["end"] <= max_span["end"]:
            keep_center.append(intv)
        if intv["start"] > max_span["end"]:
            keep_dn.append(intv)
        if intv["end"] < max_span["start"]:
            keep_up.append(intv)
        else:
            up_intv_copy = intv["start"], intv["start"] + intv["period"]
            dn_intv_copy = intv["end"] - intv["period"], intv["end"]
            # Might (should) be an unnessary first check, just needs mutual
            if intv["start"] <= max_span["end"] < intv['end'] and mutual_contain(up_intv_copy, dn_mspan_copy):
                max_span["ovl_flag"] |= 16
                intv["ovl_flag"] |= 8
                dn_need_to_choose.append(intv)
            if intv["start"] < max_span["start"] <= intv["end"] and mutual_contain(dn_intv_copy, up_mspan_copy):
                max_span["ovl_flag"] |= 8
                intv["ovl_flag"] |= 16
                up_need_to_choose.append(intv)

    if dn_need_to_choose:
        keep_dn.append(staggered_chooser(max_span, dn_need_to_choose, up=False))
    if up_need_to_choose:
        keep_up.append(staggered_chooser(max_span, up_need_to_choose)) 

    for i in simplify_region(keep_up):
        yield i
    for i in simplify_region(keep_center):
        yield i
    for i in simplify_region(keep_dn):
        yield i

def annotate_iso(annos):
    for i in annos:
        i['ovl_flag'] |= 1
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
        alt_seq = anno['motif'] * int(math.floor(anno['copies']))
        remain = int((anno['copies'] % 1) * len(anno))
        alt_seq += anno['motif'][:remain]
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
    for key in HEADER[:-1]:
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

def get_buff(reg, annos, f_start, l_end):
    """
    Try to put at least 25bp of buffer on the ends.
    """
    tree = IntervalTree()
    for anno in annos:
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

def get_interspersed(region, repeats):
    """
    Get the max scoring interspersed repeat annotation
    """
    key = f"{region['chrom']}:{region['start']}-{region['end']}"
    ret = [_[1] for _ in sorted(repeats[key]) if _[1] not in ["Simple_repeat", "Low_complexity", "Satellite", "Unknown"]]
    if ret:
        return ret[-1]
    return '.'

class AnnoTree():
    """
    Simple annotation overlap helper
    """
    def __init__(self, fn):
        self.tree, _ = truvari.build_anno_tree(fn)
        self.data = pd.read_csv(fn, sep='\t', header=None)
        self.data.columns = ['chrom', 'start', 'end', 'anno']
        def clean_anno(x):
            if x.startswith("ID"):
                return x.split(';')[0][len("ID="):]
            return x
        self.data['anno'] = self.data['anno'].apply(clean_anno)
    
    def get_annotation(self, region):
        """
        Returns the data from annotations that the region.
        Will assert that there aren't multiple overlaps.
        Returns '.' if there are no annotations
        """
        hits = self.tree[reg['chrom']].overlap(reg['start'], reg['end'])
        if not hits:
            return '.'

        if len(hits) == 1:
            return self.data.iloc[list(hits)[0].data]['anno']
        sys.stderr.write("Multianno\n%s  -  %s\n" % (reg['chrom'], str(hits)))
        ret = []
        for i in hits:
            ret.append(self.data.iloc[i.data]['anno'])
        return ",".join(ret)

def filter_homopolymers(region):
    """
    Return all of the non-homopolymer variants
    And an interval tree of where all the homopolymer runs were
    add hom_span - number of bases spanned by homopolymers
    """
    non_hom = []
    hom = IntervalTree()
    for anno in region["annos"]:
        rlen = len(anno['repeat'])
        if rlen == 1:
            hom.addi(anno['start'], anno['end'])
        else:
            non_hom.append(anno)
    hom.merge_overlaps() # Don't think I need this, but just in case
    return non_hom, hom

def resolve_homspan(region, hom_tree):
    """
    how many of the homopolymers are inside my new region bounds (but not overlapped by other annotations)
    """
    tot_bases = 0
    for i in hom_tree.overlap(region['start'], region['end']):
        if i.begin <= region['start'] and region['end'] <= i.end:
            return region['end'] - region['start']
        if i.overlaps(region['start']):
            tot_bases += i.end - region['start']
        elif i.overlaps(region['end']):
            tot_bases += region['end'] - i.begin
        else:
            tot_bases += i.length()
    return tot_bases

if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] == "0":
            for reg in iter_tr_regions("/dev/stdin"):
                print('-' * 20, reg['start'], reg['end'])
                viz_region(reg)
                print('-' * 20)
        elif sys.argv[1] == "1":
            for reg in iter_tr_regions_1("/dev/stdin"):
                print('-' * 20, reg['start'], reg['end'])
                viz_region(reg)
                print('-' * 20)
        sys.exit(0)
    in_anno, in_ref, in_rep, patho, codis = sys.argv[1:]

    reference = pysam.FastaFile(in_ref)
    repeats = joblib.load(in_rep)

    patho = AnnoTree(patho)
    codis = AnnoTree(codis)

    removed = 0
    total = 0
    for reg in iter_tr_regions(in_anno):
        total += 1
        in_anno_cnt = len(reg['annos'])

        non_hom, hom_tree = filter_homopolymers(reg)
        if not non_hom:
            removed += 1
            continue

        out_annos = simplify_region(non_hom)
        # Consolidate overlap flag
        reg['ovl_flag'] = 0
        for i in out_annos:
            reg['ovl_flag'] |= i['ovl_flag']
        
        # Do this before changing the coords because we look up by the key
        reg['interspersed'] = get_interspersed(reg, repeats)

        # Buffer work
        first_start, last_end = get_bounds(out_annos)
        new_start, new_end, s_buff, e_buff = get_buff(reg, non_hom, first_start, last_end)
        reg['start'], reg['end'] = new_start, new_end
        reg['up_buff'], reg['dn_buff'] = s_buff, e_buff
        
        # update annotations
        reg['annos'] = out_annos
        reg['n_filtered'] = in_anno_cnt - len(out_annos)
        reg['n_annos'] = len(out_annos)
        reg['n_subregions'] = len([_ for _ in get_subregions(out_annos)])
        reg['mu_purity'] = annotate_purity(reg)
        reg['pct_annotated'] = pct_annotated(reg, first_start, last_end)

        # Do it after we've update the positions
        reg['patho'] = patho.get_annotation(reg)
        reg['codis'] = codis.get_annotation(reg)

        reg["hom_span"] = resolve_homspan(reg, hom_tree)

        write_region(reg)
    sys.stderr.write(f"removed {removed} regions from {total}\n")
