"""
"""
import re
import sys
import joblib
import logging
import truvari
import pandas as pd



def do_intersection(a, b, ro=None):
    """
    Returns the intersection counts as a dataframe
    """
    cmd_tmpl = 'bedtools intersect -a {a} -b {b} -c {ro}'
    cmd_tmpl += '| cut -f4 | sort -n | uniq -c'
    if ro:
        ro = f"-r -f 0.50"
    else:
        ro = ""
    cmd = cmd_tmpl.format(a=a, b=b, ro=ro)
    ret = truvari.cmd_exe(cmd)
    data = []
    for line in ret.stdout.strip().split('\n'):
        row = re.split('\s+', line.strip())
        data.append(row)
    data = pd.DataFrame(data, columns=["count", "intersection"])
    data = data.astype(int)
    return data

if __name__ == '__main__':
    tr_bed = 'data/tr_annotated.bed'
    sources = ['baylor', 'giab', 'pacbio', 'ucsd1', 'ucsd2']
    parts = []
    truvari.setup_logging()
    for i in sources:
        logging.info("intersecting %s", i)
        dat = do_intersection(f"data/{i}/merged.bed.gz", tr_bed, ro=False)
        dat['source'] = i
        dat['ro'] = False
        parts.append(dat)
        dat = do_intersection(f"data/{i}/merged.bed.gz", tr_bed, ro=True)
        dat['source'] = i
        dat['ro'] = True
        parts.append(dat)
    data = pd.concat(parts)
    joblib.dump(data, 'data/intersection.jl')


