import sys
import gzip
import pysam
import truvari
import pandas as pd
def optional_compressed_fh(in_fn):
    """
    chooses file handler for assumed plain-text files or `*.gz` files.
    returns a generator which yields lines of a file
    """
    def gz_hdlr(fn):
        with gzip.open(fn) as fh:
            for line in fh:
                yield line.decode()

    def fh_hdlr(fn):
        with open(fn) as fh:
            for line in fh:
                yield line
    if in_fn.endswith('.gz'):
        return gz_hdlr(in_fn)
    return fh_hdlr(in_fn)


def main(in_bed, in_vcf, out_name):
    """
    """
    variants = pysam.VariantFile(in_vcf)
    with open(f"counts_{out_name}", 'w') as fout:
        for line in optional_compressed_fh(in_bed):
            chrom, start, end = line.strip().split('\t')[:3]
            start = int(start)
            end = int(end)
            cnt = 0
            bases = 0
            for i in variants.fetch(chrom, int(start), int(end)):
                # check only svs.. take this out for core analysis but keep in for extra analysis
                if 'SVLEN' not in i.info or i.info["SVLEN"] < 50:
                    continue
                vs, ve = truvari.entry_boundaries(i)
                if start <= vs and ve <= end:
                    cnt += 1
                    bases += truvari.entry_size(i)
            fout.write(f"{chrom}\t{start}\t{end}\t{cnt}\t{bases}\n")
    data = pd.read_csv(f"counts_{out_name}", sep='\t', header=None, names=['chrom', 'start', 'end', 'num_vars', 'num_bases'])
    print("statistic\tcount\tpercent")

    tot = len(data)
    print("total regions\t%d\t%d" % (tot, 1))

    no_var = data['num_vars'] == 0
    i = no_var.sum()
    print("no variant\t%d\t%.4f" % (i, i / tot))

    single_snp = (data['num_vars'] == 1) & (data['num_bases'] == 1)
    i = single_snp.sum()
    print("only a SNP\t%d\t%.4f" % (i, i / tot))

    only_snps = (data['num_vars'] > 1) & (data['num_vars'] == data['num_bases'])
    i = only_snps.sum()
    print("only SNPs\t%d\t%.4f" % (i, i / tot))

    remaining = data[~no_var & ~single_snp & ~only_snps]
    i = len(remaining)
    print("remaining\t%d\t%.4f" % (i, i / tot))

    remaining.to_csv(f"filtered_{out_name}", sep='\t', header=False, index=False)

if __name__ == '__main__':
    main(*sys.argv[1:])
