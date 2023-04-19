import pandas as pd

#strawman = pd.read_csv("../chr20_strawman_files_Feb.28.2023/strawman_regions.bed.gz",
#                       sep='\t', names=["chrom", "start", "end"])

strawman = pd.read_csv("all_refine.bed", sep="\t",
                    names=["chrom", "start", "end", "ei", "li", "th"])
strawman['key'] = strawman['ei'] + '_' + strawman['li'] + '_' + strawman['th']

states = pd.read_csv("states.txt", sep='\t')
states['key'] = states['ei'] + '_' + states['li'] + '_' + states['th']
state_map = dict(zip(states['key'], states['tier']))

strawman['tier'] = strawman['key'].map(state_map)

t1 = strawman[strawman['tier'] == 'Tier1']
t1[["chrom", "start", "end", "key"]].to_csv("Tier1.bed", sep='\t', header=False, index=False)
t2 = strawman[strawman['tier'] == 'Tier2']
t2[["chrom", "start", "end", "key"]].to_csv("Tier2.bed", sep='\t', header=False, index=False)
#print(strawman.join(states, on=["chrom", "start", "end"
