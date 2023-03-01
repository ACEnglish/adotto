import pandas as pd

d = pd.read_csv("all_merged.bed", sep='\t', names=["chrom", "start", "end", "prog"])

y = d['prog'].str.split(',').apply(lambda x: ['bench' in x, 'gangstr' in x, 'hipstr' in x, 'trgt' in x])
p_cols = ["bench", "gangstr", "hipstr", "trgt"]
y2 = pd.DataFrame(y.tolist(), index=d.index, columns=p_cols)
j = d.join(y2)

m_counts = j.groupby(p_cols).size()
print(m_counts)

from upsetplot import plot
plot(m_counts)  
from matplotlib import pyplot
pyplot.savefig('upset.png')
