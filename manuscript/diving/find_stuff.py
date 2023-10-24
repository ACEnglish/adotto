import pandas as pd

extract = (pd.read_csv("toextract.bed", sep='\t')
            .rename(columns={"chr":"chrom"})
            .set_index(["chrom", "start", "end"])
           )
fn = "/users/u233287/scratch/code/adotto/benchmark/GIABTR_benchmark.6.26/GIABTR.HG002.benchmark.regions.bed.gz"
bench = (pd.read_csv(fn, sep='\t', names=["chrom", "start", "end", "tier", 'repl', 'a', 'b', 'c','d'])
            .set_index(["chrom", "start", "end"])
            )
data = extract.join(bench)

import joblib
joblib.dump(data, 'lookat.jl')
    
