import sys
import random
import truvari

nucs = ["A", "T", "C", "G"]
def create_tandem_repeat(N=5, C=4):
    """
    Makes a random motif of length N and repeats it N times
    """
    motif = "".join(random.choices(nucs, k=N))
    sequence = motif * C
    return motif, sequence

def create_alternate(sequence, motif, copy_number, divergence=0):
    """
    create an alternate allele at a random position with some percent of bases divergent
    """
    # Extend the motif
    rep1_alt = motif * copy_number

    # the second representation will be inserted at a random position downstream
    pos = random.randint(1, len(sequence))
    rep2_alt = list(rep1_alt)

    # for the first representation, we put it at the beginning
    # now we can make the 'refence context' for w.r.t. rep2
    rep1 = rep1_alt + sequence[:pos]

    # and sometimes it'll be given some divergence
    to_change = list(range(len(rep1_alt)))
    random.shuffle(to_change)
    to_change = to_change[:int(round(len(rep1_alt) * divergence))]
    num_changes = 0
    for p in to_change:
        rep2_alt[p] = random.choice([_ for _ in nucs if _ != rep2_alt[p]])
        num_changes += 1
    rep2_alt = "".join(rep2_alt)

    # the expanded copy needs to be rotated to be placed in the new spot
    f = pos % motif_len
    rep2_alt = rep2_alt[f:] + rep2_alt[:f]
    rep2 = sequence[:pos] + rep2_alt
    return pos, rep1, rep1_alt, rep2, rep2_alt, num_changes

# Do it N times
print("motif_len motif copy_number copy_gain r1_hap r2_hap alt_pos r1_alt r2_alt div ref_context_sim unroll_sim nonroll_sim num_changes div_norm")
already_made = {} # for removing duplicates
for i in range(150):
    for motif_len in range(2, 20):
        for ref_copy_number in range(2, 50):
            motif, ref_seq = create_tandem_repeat(motif_len, ref_copy_number)
            for copy_number in range(1, 15):
                for div in range(0, 31, 5):
                    pos, rep1, rep1_alt, rep2, rep2_alt, num_changes = create_alternate(ref_seq, motif, copy_number, div / 100)
                    
                    # duplicates
                    key = rep1 + ' ' + rep2
                    if key in already_made:
                        continue
                    already_made[key] = 1

                    # The similarity using edlib is editDistance(a, b) / (len(a) + len(b))
                    # normalize the divergence of the sequence where N% of the bases in on copy is changed
                    # to account for this difference and make it a similarity
                    div_norm = 1 - (num_changes / (len(rep2_alt) * 2))
                    ref_context_sim = truvari.seqsim(rep1, rep2)
                    unroll_sim = truvari.unroll_compare(rep1_alt, rep2_alt, pos % len(rep2_alt))
                    nonroll_sim = truvari.seqsim(rep1_alt, rep2_alt)
                    # 3 TGT 2 1 TGTT TGTT 1 TGT GTT 5 1.0 0.6666666666666666 0 1.0
                    print(motif_len, motif, ref_copy_number, copy_number, 
                          rep1, rep2, pos, rep1_alt, rep2_alt, div,
                          ref_context_sim, unroll_sim, nonroll_sim, num_changes, div_norm)
