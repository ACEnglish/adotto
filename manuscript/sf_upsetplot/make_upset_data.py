import glob
import gzip
import joblib
from pathlib import Path

memberships = {}

with gzip.open("data/tr_regions.bed.gz") as fh:
    for line in fh:
        memberships[line.decode().strip()] = []

for i in glob.glob("temp/*.txt"):
    i = Path(i)
    name = i.stem
    with i.open() as fh:
        for line in fh:
            memberships[line.strip()].append(name)

joblib.dump(memberships, "memberships.jl")



    
