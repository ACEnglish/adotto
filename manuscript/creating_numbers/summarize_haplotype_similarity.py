import pandas as pd

d = pd.read_csv("haplotype_similarity.txt.gz", sep='\t', names=["chrom", "start", "end", 'a1', 'a2'])

view = pd.concat([d['a1'], d['a2']])

view = view[view != -1]

cnt = len(view)
THRESH = 0.99
over = (view < THRESH).sum()
print('total haplotypes', cnt)
print(f'haps lt {THRESH}', over, over / cnt)
mask = d['a1'].between(0, THRESH) | d['a2'].between(0, THRESH)
cnt = len(d)
print(f'regs lt {THRESH}', mask.sum(), cnt, mask.sum() / cnt)
#print(len(view))

"""
total haplotypes 177718
haps lt 0.99 8703 0.04897084144543603
regs lt 0.99 7572 102439 0.07391716045646678
"""
