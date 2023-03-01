import sys
import json
import truvari
orig_fn = sys.argv[1]
anno_fn = sys.argv[2]


feat, biotype = json.load(open(anno_fn, 'r'))

fh = truvari.opt_gz_open(orig_fn)
for line in fh:
    data = line.strip().split('\t')
    chrom, start, end = data[:3]
    key = f"{chrom}:{start}-{end}"
    if key in feat:
        f = str(feat[key])
        b = ','.join(biotype[key])
    else:
        f = "0"
        b = '.'
    data.insert(15, b)
    data.insert(15, f)
    print('\t'.join(data))
