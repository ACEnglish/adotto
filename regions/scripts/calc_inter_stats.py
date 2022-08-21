import pysam
import truvari
import itertools
import pandas as pd
from joblib import Parallel, delayed, dump
from truvari.annotations.trf import iter_tr_regions

# sys.argv[1:]
tr_regions_fn = "adotto_TRannotations_v0.2.bed.gz"
vcf_fn = "/users/u233287/scratch/code/adotto/variants/data/annotated.variants.vcf.gz"
threads = 18

def stats_counter(region_chunks):
    v = pysam.VariantFile(vcf_fn)
    ret = []
    for region in region_chunks:
        anno_count = 0
        anno_span = 0
        variant_count = 0
        trf_variant_count = 0
        known_repeat_variant_count = 0
        known_motifs = set()
        for an in region['annos']:
            anno_count += 1
            anno_span += an['end'] - an['start']
            known_motifs.add(an['repeat'])
        
        for entry in v.fetch(region['chrom'], region['start'], region['end']):
            if truvari.entry_size(entry) >= 5:
                variant_count += 1
            if 'TRFrepeat' in entry.info:
                trf_variant_count += 1
                if entry.info['TRFrepeat'] in known_motifs:
                    known_repeat_variant_count += 1
        ret.append([anno_count, anno_span, variant_count, trf_variant_count, known_repeat_variant_count])
    return ret

reg = list(iter_tr_regions(tr_regions_fn))

chunked_list = list()
chunk_size = len(reg) // threads
for i in range(0, len(reg), chunk_size):
    chunked_list.append(reg[i:i+chunk_size])

results = Parallel(n_jobs=threads, prefer="threads")(delayed(stats_counter)(i) for i in chunked_list)
m = itertools.chain(*results)
data = pd.DataFrame(m, columns=["anno_count", "anno_span", "variant_count", "trf_variant_count", "known_repeat_variant_count"])
dump(data, 'results.jl')
