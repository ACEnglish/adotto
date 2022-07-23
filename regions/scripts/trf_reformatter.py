"""
samtools faidx is 1-based

The bedfiles are (presumably) 0 based
And so are trf output is 0-based, I think.

0-based index - fetched 1-based region so I should add 1 to the starts when parsing the header

this is complex. Let's just take the header and add the start/end to the header start position to translate it
Then  I can do some pysam.FastaFile.fetch to check that the coordinates are right
"""
import sys
import joblib
import truvari
import pandas as pd


def parse_trf_output(tr_fn):
    """
    Reads the otput from tandem repeat finder
    translates the coordinates from the 'samtools faidx' fetched sequence header 
    back to whole genome coordinates
    """
    trf_cols = [("start", int),
                ("end", int),
                ("period", int),
                ("copies", float),
                ("consize", int),
                ("pctmat", int),
                ("pctindel", int),
                ("score", int),
                ("A", int),
                ("C", int),
                ("G", int),
                ("T",  int),
                ("entropy", float),
                ("repeat", str),
                ("upflank", str),
                ("sequence", str),
                ("dnflank", str)]
    with open(tr_fn, 'r') as fh:
        name = fh.readline()
        if name == "":  # no hits
            return
        name = name.strip()[1:]
        chrom, coords = name.split(':')
        wgs_start, wgs_end = coords.split('-')
        wgs_start = int(wgs_start) - 1 # 0-based correction
        wgs_end = int(wgs_end)
        while True:
            line = fh.readline()
            if line == "":
                break
            if line.startswith("@"):
                name = line.strip()[1:]
                chrom, coords = name.split(':')
                wgs_start, wgs_end = coords.split('-')
                wgs_start = int(wgs_start) - 1 # 0-based correction
                wgs_end = int(wgs_end)
                continue
            line = line.strip().split(' ')
            data = {x[0]: x[1](y) for x, y in zip(
                trf_cols, line) if x[1] is not None}
            data['chrom'] = chrom
            data['in_region_start'] = wgs_start
            data['in_region_end'] = wgs_end
            data['start'] += wgs_start
            data['end'] += wgs_start
            yield data

if __name__ == '__main__':
    #def parse_trf_output(fn):
    data = []
    out_name = sys.argv[2]
    for i in parse_trf_output(sys.argv[1]):
        data.append(i)
    #temp while I figure out how to translate this
    data = pd.DataFrame(data)
    print(truvari.optimize_df_memory(data))
    joblib.dump(data, out_name, compress=5)
