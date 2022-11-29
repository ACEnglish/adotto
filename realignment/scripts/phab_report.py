"""
Given a phab region directory, (or set of them), maybe a whole phab input.

Pull out:
How many input variants/bases
Number of variant start positions
The AFs of the input variants
Input variants mendelian error (with refhom and without)

Consistency by variants/bases of replicates


Repeat for output.vcf.gz
"""
import os
import sys
import glob
import joblib
import pysam
import truvari
import numpy as np
from collections import defaultdict, Counter

MINSIZE=5
# 3-tuple of mother, father, proband sample names
TRIOS = [("HG00513", "HG00512", "HG00514"),
         ("HG00732", "HG00731", "ei:HG00733"),
         ("NA19238", "NA19239", "ei:NA19240")
        ]

# tuple of (samp1, samp2)
REPLS = [("NA19240", "ei:NA19240"),
         ("HG00733", "ei:HG00733"),
         ("HG00733", "li:HG00733"),
         ("ei:HG00733", "li:HG00733"),
         ("HG02818", "ei:HG02818"),
         ("HG03486", "ei:HG03486"),
         ("NA12878", "li:NA12878"),
         ("NA24385", "li:NA24385")
        ]

def check_mendelian(entry, father, mother, proband):
    """
    Return 0 if family genotype tuples are mendelian consistent, 1 if not, None if there are missing genotypes
    """
    fa = entry.samples[father]["GT"]
    mo = entry.samples[mother]["GT"]
    pr = entry.samples[proband]["GT"]
    if None in fa or None in mo or None in pr:
        return None
    if pr[0] in fa and pr[1] in mo:
        return 0
    if pr[1] in fa and pr[0] in mo:
        return 0
    return 1

def update_mend(arr, entry, size, check_pass=True):
    """
    Update mendelian error counts
    """
    for idx, (m, f, p) in enumerate(TRIOS):
        is_pres = (1 in entry.samples[m]["GT"]) or (1 in entry.samples[f]["GT"]) or (1 in entry.samples[p]["GT"])
        if not is_pres:
            continue
        if check_pass:
            is_pass = entry.samples[m]["FT"] == 'P' and entry.samples[f]["FT"] == 'P' and entry.samples[p]["FT"] == 'P'
            if not is_pass:
                arr[idx][4] += 1
                continue
        arr[idx][0] += 1
        arr[idx][1] += size
        if check_mendelian(entry, f, m, p):
            arr[idx][2] += 1
            arr[idx][3] += size

def update_consis(arr, entry, size, check_pass=True):
    """
    Update consistency counts
    """
    for idx, (s1, s2) in enumerate(REPLS):
        # If present in either, update 0/1
        if not 1 in entry.samples[s1]["GT"] and not 1 in entry.samples[s2]["GT"]:
            continue
        if check_pass:
            is_pass = entry.samples[s1]["FT"] == 'P' and entry.samples[s2]["FT"] == 'P'
            if not is_pass:
                arr[idx][4] += 1
                continue
        # Present in at least one
        arr[idx][0] += 1
        arr[idx][1] += size
        if truvari.get_gt(entry.samples[s1]["GT"]) == truvari.get_gt(entry.samples[s2]["GT"]):
            arr[idx][2] += 1
            arr[idx][3] += size

def var_filter(entry, reg_st, reg_ed):
    """
    Is this an entry to analyze
    """
    st, ed = truvari.entry_boundaries(entry)
    if not (st >= reg_st and ed <= reg_ed):
        return True, None
    m_size = truvari.entry_size(entry)
    if m_size < MINSIZE:
        return True, m_size
    return False, m_size

def split_region(fn):
    """
    split chr:start-end to their pieces
    """
    name = os.path.dirname(fn).split('/')[-1]
    c, coord = name.split(':')
    reg_st, reg_ed = [int(x) for x in coord.split('-')]
    return c, reg_st, reg_ed

def vcf_stats(fn, check_pass=True):
    """
    Get stats from a VCF.
    """
    if not os.path.exists(fn):
        return "FileNotFound"
    var_count = 0
    var_bases = 0
    var_positions = set()
    allele_counts = []
    # Need to remove variants that are part of 'base'
    # only because of the initial fetch but don't have start/end within region

    _, reg_st, reg_ed = split_region(fn)

    # Rows are each trio
    # count/bases in trio, c/b for all errors
    menderr = np.zeros((len(TRIOS), 5))
    # Rows are each pair
    # count/bases for nonref and consistent
    consis = np.zeros((len(REPLS), 5))
    
    vcf = pysam.VariantFile(fn)
    for entry in vcf:
        flt, m_size = var_filter(entry, reg_st, reg_ed)
        if flt:
            continue
        var_count += 1
        if m_size == 0: # SNPs and MNPs
            m_size = len(entry.ref)
        var_bases += m_size
        var_positions.add(entry.start)
        allele_counts.append(truvari.calc_af([_["GT"] for _ in entry.samples.values()])["AC"][1])
        update_mend(menderr, entry, m_size, check_pass)
        update_consis(consis, entry, m_size, check_pass)
    
    var_positions = len(var_positions)
    return {"var_count": var_count,
            "var_bases": var_bases,
            "var_positions": var_positions,
            "allele_counts": allele_counts,
            "menderr": menderr,
            "consis": consis}
 
def tr_stats(fn):
    """
    Collect the motif counts
    Should be 1-to-1 with AC if that makes a difference
    """
    motif_counts = [] # Motif TRFdiff (may be none)
    if not os.path.exists(fn):
        return "FileNotFound"
    _, reg_st, reg_ed = split_region(fn)
    try:
        fh = pysam.VariantFile(fn)
    except ValueError:
        return "FileCorrupted"
    for entry in fh:
        flt, m_size = var_filter(entry, reg_st, reg_ed)
        if flt:
            continue
        if "TRFrepeat" in entry.info:
            rep = entry.info['TRFrepeat']
        else:
            rep = None
        if "TRFdiff" in entry.info:
            diff = round(entry.info["TRFdiff"], 1) # floats returned by pysam sometimes get wonky
        else:
            diff = None
        motif_counts.append([rep, diff])
    return motif_counts

def full_dir():
    in_dir = sys.argv[1]
    results = {}
    for region_dir in glob.glob(os.path.join(in_dir, "*:*-*")):
        #for region in glob.glob(in_dir):
        m_name = os.path.basename(region_dir.rstrip('/'))
        results[m_name] = {'in': vcf_stats(os.path.join(region_dir, 'base.vcf.gz')),
                       'out': vcf_stats(os.path.join(region_dir, 'output.vcf.gz'), check_pass=False),
                       'inTR': tr_stats(os.path.join(region_dir, 'base.vcf.gz')),
                       'outTR': tr_stats(os.path.join(region_dir, 'output.trf.vcf'))}

    joblib.dump(results, f'phabsummary_{MINSIZE}bp.jl')

def single_dir():
    region_dir = sys.argv[1]
    results = {}
    m_name = os.path.basename(region_dir.rstrip('/'))
    results[m_name] = {'in': vcf_stats(os.path.join(region_dir, 'base.vcf.gz')),
                       'out': vcf_stats(os.path.join(region_dir, 'output.vcf.gz'), check_pass=False),
                       'inTR': tr_stats(os.path.join(region_dir, 'base.vcf.gz')),
                       'outTR': tr_stats(os.path.join(region_dir, 'output.trf.vcf'))}
    print(results)

if __name__ == '__main__':
    #single_dir()
    full_dir()

