import pandas as pd

d = pd.read_csv("Finished_v1_bench_regions.bed", sep='\t', 
                names=['chrom', 'start', 'end', 'tier', 'repl', 'vflag', 'entropy', 'ad1', 'ad2'])
d['entropy'] = d['entropy'].round(4)
d.to_csv("Finished_v1_bench_regions_round_entropy.bed", sep='\t', header=False, index=False)
