import sys
import pandas as pd

d = pd.read_csv(sys.argv[1], sep='\t')


#cols = ['HG002', "NA24385", "li:NA24385"]
cols = [_.strip() for _ in open('males.txt', 'r')] # males
looking_for = ['1,1,0,0', '0,0,1,1']
#v = (~d[cols].isin(looking_for)).sum()

#cols = d.columns[4:]
#looking_for = ['1,1,1,1']
v = (~d[cols].isin(looking_for)).sum()
# sample, inadequate counts
print(v.value_counts())

print(v)
idx = v[v != 0]
print(len(idx), len(v))
#print(idx)
