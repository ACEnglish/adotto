import glob

files = glob.glob("logs/*")
rows = []
for f in files:
    if "CHM13Y" in f: continue

    if f.startswith("logs/hprc"):
        d = f.split('.')[1:-2]
        d.insert(1, 'hprc')
        with open(f) as fh:
            d.extend(fh.readline().strip().split(' '))
            rows.append(d)
    else:
        with open(f) as fh:
            rows.append(fh.readline().strip().split(' '))
import pandas as pd
d = pd.DataFrame(rows, columns=["sample", "project", "haplotype", "n_contigs", "n50"])
d.to_csv("assembly_stats.txt", sep='\t', index=False)
        

            
        
