import sys
import pysam
import truvari
import pandas as pd
import joblib
import logging
from collections import defaultdict
truvari.setup_logging()
vcf_fn=sys.argv[1]
tree_fn = sys.argv[2]
table_fn = sys.argv[3]

tree = joblib.load(tree_fn)
cov = joblib.load(table_fn)
cov.sort_index(inplace=True)

vcf = pysam.VariantFile(vcf_fn)
n_header = vcf.header.copy()
#n_header.add_line(('##FORMAT=<ID=BPDP,Number=.,Type=Integer,'
#                   'Description="Coverage of Assemblies of start/end breakpoints">'))
#n_header.add_line(('##FORMAT=<ID=AD,Number=R,Type=Integer,'
#                   'Description="Coverage of variant (from BPDP)">'))
#n_header.add_line(('##FORMAT=<ID=FT,Number=1,Type=String,'
#                   'Description="Genotype passes for having diploid single coverage or fails">'))


out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
for entry in vcf:
    up_pos_coverage = cov.loc[[_ for _ in list(tree[entry.chrom][entry.pos])[0].data]]['cov']
    up_pos_coverage.index = [".".join(_.split('.')[:2]) for _ in up_pos_coverage.index]
    dn_pos_coverage = cov.loc[[_ for _ in list(tree[entry.chrom][entry.pos])[0].data]]['cov']
    dn_pos_coverage.index = [".".join(_.split('.')[:2]) for _ in dn_pos_coverage.index]

    for sample in entry.samples:
        if sample + '.0' not in up_pos_coverage:
            logging.warning("missing up %s.0 at %s:%d.0", sample, entry.chrom, entry.start)
        if sample + '.1' not in up_pos_coverage:
            logging.warning("missing up %s at %s:%d.1", sample, entry.chrom, entry.start)
        if sample + '.0' not in dn_pos_coverage:
            logging.warning("missing dn %s at %s:%d.0", sample, entry.chrom, entry.start)
        if sample + '.1' not in dn_pos_coverage:
            logging.warning("missing dn %s at %s:%d.1", sample, entry.chrom, entry.start)

        if entry.samples[sample]["FT"] != '.':
            # Trust the earlier step got the coverage right
            continue
        
        u_cov1 = int(up_pos_coverage[sample + '.0'])
        u_cov2 = int(up_pos_coverage[sample + '.1'])
        d_cov1 = int(dn_pos_coverage[sample + '.0'])
        d_cov2 = int(dn_pos_coverage[sample + '.1'])
 
        entry.samples[sample]["BPDP"] = [u_cov1, d_cov1, u_cov2, d_cov2]
        is_single1 = u_cov1 == 1 and d_cov1 == 1 
        is_single2 = u_cov2 == 1 and d_cov2 == 1
        if is_single1 and is_single2:
            entry.samples[sample]["FT"] = 'PASS'
        else:
            entry.samples[sample]["FT"] = 'FAIL'
        # Set the genotype
        gt1 = 0 if u_cov1 >= 1 and d_cov1 >= 1 else None
        gt2 = 0 if u_cov2 >= 1 and u_cov2 >= 1 else None
        entry.samples[sample]["GT"] = (gt1, gt2)
        entry.samples[sample]["AD"] = (max(u_cov1, d_cov1), max(u_cov2, d_cov2))
    out.write(entry)

logging.info("Finished annotate_dip_cov")
