import truvari
import json
import sys

data = json.load(open(sys.argv[1]))
for i in ["DEL", "INS"]:
    for j in truvari.SZBINS:
        if j in data[i]:
            print(i, j, data[i][j], sep='\t')
