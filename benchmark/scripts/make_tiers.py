import sys
import pysam
import truvari
import numpy as np
import pandas as pd
from truvari.annotations.lcr import sequence_entropy

def get_var_state(vcf, chrom, start, end):
    """
    Return the variant state of region
    fivebp_state - HG002, Other, Both (note- we gotta ignore the HG002 replicates)
    small_state - HG002, Other, Both (note- we gotta ignore the HG002 replicates)
    0 - No variant
    1 - HG002 >=5bp variant
    2 - HG002 <5bp variant
    4 - Other >=5bp variant
    8 - Other <5bp variant
    # so max'd out is 15
    """
    presence = 0
    for entry in vcf.fetch(chrom, start, end):
        # Stop when we've seen everything
        if presence == 15:
            break
        st, ed = truvari.entry_boundaries(entry)
        # must be within to count
        if not (start <= st <= ed <= end):
            continue
        is_hg002 =  1 in entry.samples["HG002"]["GT"]
        is_other = False
        for samp in entry.samples:
            if is_other: #No more reason to check
                continue
            if samp in ["HG002", "li:NA24385", "NA24385"]:
                continue
            if 1 in entry.samples[samp]["GT"]:
                is_other = True
        sz = truvari.entry_size(entry)
        if sz >= 5:
            presence |= 1 if is_hg002 else 0
            presence |= 4 if is_other else 0
        else:
            presence |= 2 if is_hg002 else 0
            presence |= 8 if is_other else 0
    return presence

def main():
    refine_bed = sys.argv[1]
    curated_states = sys.argv[2]
    ref_fn = sys.argv[3]
    vcf_fn = sys.argv[4]
    out_fn = sys.argv[5]

    strawman = pd.read_csv(refine_bed, sep="\t",
                        names=["chrom", "start", "end", "ei", "li", "th"])
    strawman['reps'] = strawman['ei'] + '_' + strawman['li'] + '_' + strawman['th']

    states = pd.read_csv(curated_states, sep='\t')
    states['key'] = states['ei'] + '_' + states['li'] + '_' + states['th']
    state_map = dict(zip(states['key'], states['tier']))

    strawman['tier'] = strawman['reps'].map(state_map)

    ref = pysam.FastaFile(ref_fn)
    vcf = pysam.VariantFile(vcf_fn)
    all_states = []
    all_entropy = []
    #print("\t".join(["chrom", "start", "end", "replicate", "var_state", "entropy"]))
    for idx, row in strawman.iterrows():
        all_states.append(get_var_state(vcf, row["chrom"], row["start"], row["end"]))
        all_entropy.append(sequence_entropy(ref.fetch(row["chrom"], row["start"], row["end"])))
    strawman["var_state"] = all_states
    strawman["entropy"] = all_entropy
    # chrom start end replicate var_state entropy
    strawman[['chrom', 'start', 'end', 'tier', 'reps', 'var_state', 'entropy']].to_csv(out_fn, sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()
