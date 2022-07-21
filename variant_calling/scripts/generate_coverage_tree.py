"""
Given an input file of sample\thap1_cov.bed\thap2_cov.bed
Make a giant tree for lookup
"""
import sys
import joblib
import truvari
import logging
import itertools
import pandas as pd

def parse_bed_file(fn, sample, allele):
    idxfmt = sample + '.' + allele + '.' + '{}'
    tree, cnt = truvari.build_anno_tree(fn, idxfmt=idxfmt)
    cov_info = pd.read_csv(fn, sep='\t', names=['chrom', 'start', 'end', 'cov'])
    cov_info["sample"] = sample
    cov_info["allele"] = allele
    cov_info['idx'] = cov_info.index
    def name_joiner(x):
        return f"{x['sample']}.{x['allele']}.{x['idx']}"
    cov_info['key'] = cov_info.apply(name_joiner, axis=1)
    cov_info.set_index('key', inplace=True)
    return tree, cov_info

if __name__ == '__main__':
    """
    """
    cov_files = sys.argv[1]
    truvari.setup_logging()
    cov_parts = []
    tree_parts = []
    with open(cov_files) as fh:
        for line in fh:
            sample, cov1, cov2 = line.strip().split('\t')

            logging.info("parsing %s.%s", sample, cov1)
            tree1, cov1 = parse_bed_file(cov1, sample, '0')
            cov_parts.append(cov1)
            tree_parts.append(tree1)

            logging.info("parsing %s.%s", sample, cov2)
            tree2, cov2 = parse_bed_file(cov2, sample, '1')
            cov_parts.append(cov2)
            tree_parts.append(tree2)

    logging.info("dumping coverage")
    coverage = pd.concat(cov_parts)
    joblib.dump(coverage, "coverage.jl")

    logging.info("dumping pre_trees")
    joblib.dump(tree_parts, 'pre_trees.jl')

    logging.info("tree union")
    all_chroms = set(itertools.chain(*[_.keys() for _ in tree_parts]))
    logging.info(all_chroms)
    tree = tree_parts[0]
    for i in tree_parts[1:]:
        for chrom in all_chroms:
            tree[chrom].update(i[chrom])
    logging.info("tree split/merge")
    def m_reducer(a, b):
        a.append(b)
        return a
    for chrom in tree:
        logging.info(chrom)
        tree[chrom].split_overlaps()
        tree[chrom].merge_overlaps(data_reducer=m_reducer, data_initializer=[])
        
    logging.info("dumping trees")
    joblib.dump(tree, 'annotree.jl')

    logging.info("finished")
