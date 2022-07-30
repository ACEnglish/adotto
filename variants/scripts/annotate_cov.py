import sys
import pysam
import truvari
import pandas as pd

vcf_fn=sys.argv[1]
bed_fn=sys.argv[2]
tree, cnt = truvari.build_anno_tree(bed_fn)

cov_info = pd.read_csv(bed_fn, sep='\t', names=['chrom', 'start', 'end', 'cov'])
vcf = pysam.VariantFile(vcf_fn)
n_header = vcf.header.copy()
n_header.add_line(('##FORMAT=<ID=BPDP,Number=.,Type=Integer,'
                   'Description="Coverage of Assemblies of start/end breakpoints">'))
n_header.add_line(('##FILTER=<ID=PASS,Description="All Filters Passed">'))
n_header.add_line(('##FILTER=<ID=COV,Number=0,Type=Flag,'
                   'Description="Coverage of variant != 1">'))

out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
hap_cnt = len(vcf.header.samples) * 2
for entry in vcf:
    entry.translate(n_header)
    up_pos = int(cov_info.iloc[list(tree[entry.chrom].at(entry.start))[0].data]['cov'])
    dn_pos = int(cov_info.iloc[list(tree[entry.chrom].at(entry.stop))[0].data]['cov'])
    entry.samples[0]["BPDP"] = [up_pos, dn_pos]
    is_single = up_pos == 1 and dn_pos == 1
    if is_single:
        entry.filter.add('PASS')
    else:
        entry.filter.add('COV')
    out.write(entry)
