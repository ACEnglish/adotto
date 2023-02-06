import sys
import truvari
remove_keys = {}
fh = truvari.opt_gz_open("blue.bed.gz")
for line in fh:
    line = line.strip().split('\t')
    remove_keys[line[0] + ':' + line[1]] = 1
fh = truvari.opt_gz_open("green.bed.gz")
for line in fh:
    line = line.strip().split('\t')
    remove_keys[line[0] + ':' + line[1]] = 1
with open("../all_chr20_covered.bed", 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')
        key = data[0] + ':' + data[1]
        if key not in remove_keys:
            sys.stdout.write(line)
