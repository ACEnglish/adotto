"""
Puts in the coverage information as well as consolidates SAMPLE[1] into SAMPLE[0]
such that cut -f1-10 on the VCF will then have a single diploid sample's GT
"""
import sys
import pysam
import truvari
import pandas as pd

vcf_fn=sys.argv[1]

cov_files=sys.argv[2]
#sample_name\tallele1\tallele2

# Load all the bed files into one big tree.

bed_fn1=sys.argv[2]
bed_fn2=sys.argv[3]
tree1, cnt1 = truvari.build_anno_tree(bed_fn1)
cov_info1 = pd.read_csv(bed_fn1, sep='\t', names=['chrom', 'start', 'end', 'cov'])
tree2, cnt2 = truvari.build_anno_tree(bed_fn2)
cov_info2 = pd.read_csv(bed_fn2, sep='\t', names=['chrom', 'start', 'end', 'cov'])


vcf = pysam.VariantFile(vcf_fn)
n_header = vcf.header.copy()
n_header.add_line(('##FORMAT=<ID=BPDP,Number=.,Type=Integer,'
                   'Description="Coverage of Assemblies of start/end breakpoints">'))
n_header.add_line(('##FORMAT=<ID=AD,Number=R,Type=Integer,'
                   'Description="Coverage of variant (from BPDP)">'))
n_header.add_line(('##FORMAT=<ID=FT,Number=1,Type=String,'
                   'Description="Genotype passes for having diploid single coverage or fails">'))


out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
hap_cnt = len(vcf.header.samples) * 2
for entry in vcf:
    entry.translate(n_header)
    up_pos1 = int(cov_info1.iloc[list(tree1[entry.chrom].at(entry.start))[0].data]['cov'])
    dn_pos1 = int(cov_info1.iloc[list(tree1[entry.chrom].at(entry.stop))[0].data]['cov'])
    up_pos2 = int(cov_info2.iloc[list(tree2[entry.chrom].at(entry.start))[0].data]['cov'])
    dn_pos2 = int(cov_info2.iloc[list(tree2[entry.chrom].at(entry.stop))[0].data]['cov'])
    entry.samples[0]["BPDP"] = [up_pos1, dn_pos1, up_pos2, dn_pos2]
    is_single1 = up_pos1 == 1 and dn_pos1 == 1 
    is_single2 = up_pos2 == 1 and dn_pos2 == 1
    entry.filter.clear()
    entry.filter.add("PASS")
    if is_single1 and is_single2:
        entry.samples[0]["FT"] = 'PASS'
    else:
        entry.samples[0]["FT"] = 'FAIL'
    #Combine the genotypes
    gt1 = entry.samples[0]["GT"][0]
    if gt1 is None and (up_pos1 >= 1 or dn_pos1 >= 1):
        gt1 = 0
    gt2 = entry.samples[1]["GT"][0]
    if gt2 is None and (up_pos2 >= 1 or dn_pos2 >= 1):
        gt2 = 0
    entry.samples[0]["GT"] = (gt1, gt2)
    entry.samples[0]["AD"] = (max(up_pos1, up_pos1), max(up_pos2, dn_pos2))
    out.write(entry)
