import sys
import gzip
import pysam
import truvari

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
    

def main(in_bed, in_vcf):
    """
    """
    variants = pysam.VariantFile(in_vcf)
    for line in optional_compressed_fh(in_bed):
        chrom, start, end = line.strip().split('\t')[:3]
        cnt = 0
        bases = 0
        for i in variants.fetch(chrom, int(start), int(end)):
            cnt += 1
            bases += truvari.entry_size(i)
        print(f"{chrom}\t{start}\t{end}\t{cnt}\t{bases}")

if __name__ == '__main__':
    main(*sys.argv[1:])
