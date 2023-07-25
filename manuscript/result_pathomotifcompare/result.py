catalog_fn = 
patho_fn = 


def smallest_roll(seq):
    if isinstance(seq, float):
        return "GCN" # HOXA13 and ARX cheat
    if ',' in seq:
        return seq
    sm = sorted(list(seq))
    i = seq.index(sm[0])
    return seq[i:] + seq[:i]

for name, i in [('A', cluster_A), ('B', cluster_B), ('C', cluster_C), ('D', cluster_D), ('E', cluster_E)]:
    x = view[i]['motif'].apply(smallest_roll).value_counts()
    display(x)
    print('-' * 10)
