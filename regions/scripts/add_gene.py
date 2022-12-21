
import sys
import json
from collections import Counter, defaultdict

class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)

def get_biotype(text):
    for i in text.split(';'):
        i = i.strip()
        if i.startswith("gene_biotype"):
            return i.split('"')[1]
    return None

feature_flags = {
         "transcript": 1,
         "gene": 1,
         "five_prime_utr": 2,
         "start_codon": 4,
         "CDS": 8,
         "exon": 16,
         "stop_codon": 32,
         "three_prime_utr": 64,
         "Selenocysteine": 128,
        }

feat_lookup = Counter()
biotype_lookup = defaultdict(set)
#cnt = 40
with open(sys.argv[1], 'r') as fh:
    for line in fh:
        #cnt -= 1
        #if cnt <= 0:
            #break
        data = line.strip().split('\t')
        if data[18] == '.':
            continue
        chrom, start, end = data[:3]
        key = f"{chrom}:{start}-{end}"
        feature = data[18]
        feat_lookup[key] |= feature_flags[feature]
        bt = get_biotype(data[24])
        if bt:
            biotype_lookup[key].add(bt)

print(json.dumps([feat_lookup, biotype_lookup], indent=4, cls=SetEncoder))
