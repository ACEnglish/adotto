import sys
import pandas as pd

d = pd.read_csv(sys.argv[1], sep='\t')

males = [_.strip() for _ in open('males.txt', 'r')]

looking_for = ['1,1,0,0', '0,0,1,1']
v = (~d[males].isin(looking_for)).sum()
# sample, inadequate counts
print(v.value_counts())

print(v)
idx = v[v != 0]
print(len(idx), len(v))
#print(idx)
