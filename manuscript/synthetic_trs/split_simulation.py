"""
Similar to seq similarity, we're going to make some random TRs.

But, for the rep2, we're going to split them up.
For example, +2 copies, we put them in two different positions.
(Do I we want to allow partials splits?)

Then we want to run them through mafft and see if variants are starting/ending
in the same position, having the same counts.

We can again do all this in-memory.

Ho - harmonized variants are equally comparable us unharmonized vriants
Ha - harmonization improves comprability

I'm also going to do it with pywfa. Its faster (and available), but highlights the need for MSA over independent
realignment.
"""
import sys
import json
import random
import pysam
import truvari
import pandas as pd
from io import BytesIO
from pywfa.align import WavefrontAligner
from truvari.phab import fasta_reader, expand_cigar
from truvari.msatovcf import msa_to_vars

nucs = ["A", "T", "C", "G"]
def substitute(seq, div):
    """
    Randomly substitute some percent of bases in a sequence
    """
    m_len = len(seq)
    num_changes = max(1, int(round(m_len * div)))
    to_change = list(range(m_len))
    random.shuffle(to_change)
    new_seq = list(seq)
    for p in to_change[:num_changes]:
        new_seq[p] = random.choice([i for i in nucs if i != new_seq[p]])
    return "".join(new_seq), num_changes

def split(seq, motif_len, num):
    """
    Randomly chop up a sequence
    """
    # the start positions of the motifs (easier rolling)
    split_motif = []
    for i in range(0, len(seq), motif_len):
        split_motif.append(seq[i:i + motif_len])
    # to make num parts we do num-1 splits
    # do the 1 to len-1 because we want at least one motif in each split
    pos = sorted(random.sample(range(1, len(seq) + 1, motif_len), k=num))
    start = 0
    for i in pos:
        yield "".join(split_motif[start:i])
        start = i
    yield "".join(split_motif[start:])

def split_insert(alt, motif_len, seq, num_splits):
    """
    Randomly split the motifs and then insert the alternate
    """
    parts = split(alt, motif_len, num_splits)
    
    # Only inserting at the beginning of one of the motifs in the squence
    m_chunks = range(motif_len, len(seq), motif_len)
    insert_at = sorted(random.sample(m_chunks, k=num_splits), reverse=True)
    n_seq = list(seq)
    for subseq, pos in zip(parts, insert_at):
        # I don't need to roll here because I'm putting them between copies
        n_seq.insert(pos, subseq)
    return "".join(n_seq)

def create_alternate(sequence, repeat, copy_number, divergence=0, n_splits=1):
    """
    We'll have a trf annotation to work with
    We'll do a simple insertion at the beginning for one
    Do the divergence stuff
    Then do the insertion stuff.
    We'll return how many split (from insertion stuff)
    And the two haplotypes
    Then send ref/two haps to mafft/pywfa, count how many variants in rep2 after
    """
    alt = repeat * copy_number
    rep1 = alt + sequence

    rep2_alt, num_changes = substitute(alt, divergence)
    rep2 = split_insert(rep2_alt, len(repeat), sequence, n_splits)

    return rep1, rep2, num_changes

def create_alternate_same(sequence, repeat, copy_number, divergence=0, n_splits=1):
    """
    We'll have a trf annotation to work with
    We'll do a simple insertion at the beginning for one
    Do the divergence stuff
    Then do the insertion stuff.
    We'll return how many split (from insertion stuff)
    And the two haplotypes
    Then send ref/two haps to mafft/pywfa, count how many variants in rep2 after
    """
    alt = repeat * copy_number
    sub_alt, num_changes = substitute(alt, divergence)
    rep1 = sub_alt + sequence
    rep2 = split_insert(rep2_alt, len(repeat), sequence, n_splits)

    return rep1, rep2, num_changes

def harmonize_wfa(ref, rep1, rep2):
    """
    just mafft or wfa of the three sequences and then counting the number of variants
    returns the set of positions where there's a variant

    Wrong. I need to expand cigar and make variants just like with mafft
    """
    aligner = WavefrontAligner(ref, span="end-to-end", heuristic="adaptive")
    # Turn this into an msa
    msa = {'ref_': ref}
    aligner.wavefront_align(rep1)
    msa['rep_1_'] = expand_cigar(rep1, ref, aligner.cigartuples)
    aligner.wavefront_align(rep2)
    msa['rep_2_'] = expand_cigar(rep2, ref, aligner.cigartuples)
    return msa_var_counter(msa)


def harmonize_mafft(ref, rep1, rep2):
    """
    just mafft or wfa of the three sequences and then counting the number of variants
    returns the set of positions where there's a variant
    """
    cmd = f"mafft --quiet --auto --thread 1 -"
    seq_bytes = BytesIO()
    seq_bytes.write(f">ref_\n{ref}\n>rep_1_\n{rep1}\n>rep_2_\n{rep2}".encode())
    seq_bytes.seek(0)
    ret = truvari.cmd_exe(cmd, stdin=seq_bytes.read())
    msa = {}
    for name, sequence in fasta_reader(ret.stdout, name_entries=False):
        msa[name] = sequence.decode().upper()
    return msa_var_counter(msa)

def msa_var_counter(msa):
    rep1_vars = set()
    rep2_vars = set()
    
    _, variants = msa_to_vars(msa, 'none', msa['ref_'], 0, 'N')
    for v, s in variants.items():
        d = v.split('\t')
        if len(d[3]) == 1 and len(d[4]) == 1: # Skip SNPs
            continue
        pos = d[1]
        if 'rep_1' in s:
            rep1_vars.add(pos)
        if 'rep_2' in s:
            rep2_vars.add(pos)
    return len(rep1_vars), len(rep2_vars), len(rep1_vars.union(rep2_vars))

def test():
    motif = "ATCG"
    ref_seq = motif * 4
    num_splits = 2
    num_copies = 2
    #motif = motif.lower()
    alt = motif * num_copies
    #rep1, rep2, nc = split_insert(alt, len(motif), ref_seq, num_splits)
    rep1, rep2, nc = create_alternate(ref_seq, motif, num_copies, 0.05, num_splits)
    v1, v2, num_pos = harmonize_mafft(ref_seq, rep1, rep2)
    print(rep1, rep2, len(motif), len(alt), nc, num_splits, v1, v2, num_pos)
    v1, v2, num_pos = harmonize_wfa(ref_seq, rep1, rep2)
    print(rep1, rep2, len(motif), len(alt), nc, num_splits, v1, v2, num_pos)

def generate_random_spots(reference, catalog, num_iter):
    """
    Over num_iter spots
    pull the reference sequence from one of the annotations
    generate between 2 and 10 copies of the motif to be inserted
    generate 0% to 30% divergence
    generate between 1 and 4 splits

    Filters:
    - no interspersed
    - high purity 95+
    """
    view = catalog[(catalog['interspersed'] == '.') & (catalog['mu_purity'] >= 95)]
    for i in range(num_iter):
        annos = json.loads(catalog.sample(1).iloc[0]['annos'])
        # get the longest spanning motif in case of multiple
        max_span = 0
        best = None
        for an in annos:
            span = an['end'] - an['start']
            if span > max_span:
                max_span = span
                best = an
        ref_seq = reference.fetch(best['chrom'], best['start'], best['end'])
        motif = best['motif']
        for num_copies in range(2, 15):
            for divergence in [0, 0.05, 0.10]:#, 0.15, 0.20, 0.25, 0.30]:
                for num_splits in range(1, 5):
                    yield ref_seq, motif, num_copies, divergence, num_splits

def main():
    ref_fn = "/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa"
    cat_fn = "/Users/english/code/adotto/regions/adotto_TRregions_v1.1.bed"
    N = 1000 #Number of loci we'll be simulating
    ref = pysam.FastaFile(ref_fn)
    cat = pd.read_csv(cat_fn, sep='\t')
    print('num_copies num_changes motif num_splits m_num_vars1 m_num_vars2 m_num_pos w_num_vars1 w_num_vars2 w_num_pos rep1 rep2')
    for reference, motif, num_copies, divergence, num_splits in generate_random_spots(ref, cat, N):
        try:
            rep1, rep2, num_changes = create_alternate(reference, motif, num_copies, divergence, num_splits)
            w_v1, w_v2, w_num_pos = harmonize_wfa(reference, rep1, rep2)
            m_v1, m_v2, m_num_pos = harmonize_mafft(reference, rep1, rep2)
            print(num_copies, num_changes, motif, num_splits, m_v1, m_v2, m_num_pos, w_v1, w_v2, w_num_pos, rep1, rep2)
        except Exception:
            # Sometimes the parameters aren't adequate for doing a split. 
            # Instead of handling every edge case, we'll just skip them.
            pass


if __name__ == '__main__':
    #test()
    main()
