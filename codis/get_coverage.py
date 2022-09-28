import sys
import pysam
import truvari
import pandas as pd
import joblib
import logging
from collections import defaultdict
truvari.setup_logging()

bed_fn=sys.argv[1]
tree_fn = sys.argv[2]
table_fn = sys.argv[3]

tree = joblib.load(tree_fn)
cov = joblib.load(table_fn)
cov.sort_index(inplace=True)

all_samples = cov['sample'].unique()
vcf = open(bed_fn, 'r')

out = open('/dev/stdout', 'w')
out.write("chrom\tstart\tend\tname\t" + "\t".join(all_samples) + "\n")
for entry in vcf:
    entry = entry.strip()
    chrom, start, end, name = entry.split('\t')
    start = int(start)
    end = int(end)
    up_pos_coverage = cov.loc[[_ for _ in list(tree[chrom][start])[0].data]]['cov']
    up_pos_coverage.index = [".".join(_.split('.')[:2]) for _ in up_pos_coverage.index]
    dn_pos_coverage = cov.loc[[_ for _ in list(tree[chrom][end])[0].data]]['cov']
    dn_pos_coverage.index = [".".join(_.split('.')[:2]) for _ in dn_pos_coverage.index]
    
    out.write(entry)
    for sample in all_samples:
        if sample + '.0' not in up_pos_coverage:
            logging.warning("missing up %s.0 at %s:%d.0", sample, entry.chrom, entry.start)
        if sample + '.1' not in up_pos_coverage:
            logging.warning("missing up %s at %s:%d.1", sample, entry.chrom, entry.start)
        if sample + '.0' not in dn_pos_coverage:
            logging.warning("missing dn %s at %s:%d.0", sample, entry.chrom, entry.start)
        if sample + '.1' not in dn_pos_coverage:
            logging.warning("missing dn %s at %s:%d.1", sample, entry.chrom, entry.start)

        u_cov1 = int(up_pos_coverage[sample + '.0'])
        u_cov2 = int(up_pos_coverage[sample + '.1'])
        d_cov1 = int(dn_pos_coverage[sample + '.0'])
        d_cov2 = int(dn_pos_coverage[sample + '.1'])
        out.write(f"\t{u_cov1},{d_cov1},{u_cov2},{d_cov2}")
        
    out.write('\n')

logging.info("Finished annotate_dip_cov")
