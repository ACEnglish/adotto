import pandas as pd


# bedtools intersect -wao -a <(cut -f1-3 adotto_TRregions_v1.1.bed) -b <(cut -f1-4../pathogenic/Patho.tsv) > inter.bed
# Then manually curate it for sites that hit too many spots (also header)
intersect = pd.read_csv("inter.bed", sep='\t').set_index(['chr', 'start', 'end'])
catalog = pd.read_csv("adotto_TRregions_v1.1.bed", sep='\t').set_index(['chr', 'start', 'end'])

catalog['patho'] = intersect['locus']
catalog['patho'] = catalog['patho'].fillna('.')
catalog.reset_index().to_csv("new.bed", sep='\t', index=False, header=False)
# Then, make index on catalog
# Then update to the Locus for the catalog

# pyranges?
