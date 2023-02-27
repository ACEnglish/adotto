import pandas as pd

d = pd.read_csv("grch38.gene.raw.txt.gz", sep='\t')

d[~d['proteinID'].isna()][["chrom", "txStart", "txEnd"]].to_csv("grch38.proteinGenes.bed", sep='\t', header=False, index=False)
