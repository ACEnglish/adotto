import sys
import json
import tabix
import truvari

def main(reg_fn, anno_fn):
    """
    Combine regions and annotations into a file
    """
    header = [('chrom', str), ('start', int), ('end', int), ('period', float),
              ('copies', float), ('score', int), ('entropy', float), ('repeat', str)]

    tb = tabix.open(anno_fn)
    for line in truvari.opt_gz_open(reg_fn):
        chrom, start, end = line.strip().split('\t')[:3]
        start = int(start)
        end = int(end)
        m_data = []
        for i in tb.query(chrom, start, end):
            try:
                m_data.append({fmt[0]: fmt[1](x) for x, fmt in zip(i, header)})
            except tabix.TabixError as e:
                pass # allow empty regions
        data_str = json.dumps(m_data)
        print(f"{chrom}\t{start}\t{end}\t{data_str}")

if __name__ == '__main__':
    reg, anno = sys.argv[1:]
    main(reg, anno)
