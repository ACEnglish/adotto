from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta
from truvari.phab import fasta_reader
import pandas as pd
import joblib
import json
import re
import pysam

def dedup(haps, refseq):
    haps = dict(fasta_reader(haps.decode(), False))
    for k in haps:
        haps[k] = haps[k].decode().strip() #annoying
    ref_len = len(refseq)
    deduped = {} # length : sequence
    acs = {} # length: allele count
    hg002 = [None, None] # any length that hg002 has
    for hap in haps:
        # Skip replicates
        if hap.startswith("NA24385") or hap.split('_')[0].count(':'):
            continue

        m_length = len(haps[hap])
        if m_length in acs:
            acs[m_length] += 1
        else:
            acs[m_length] = 1

        if hap.startswith("HG002_"):
            if hap.startswith("HG002_1"):
                hg002[0] = m_length
            else:
                hg002[1] = m_length

        if m_length not in deduped:
            deduped[m_length] = haps[hap]
    ret = ""
    wrote_ref = False
    for length in deduped:
        if length == hg002[0] and length == hg002[1]:
            label = "HG002:"
        elif length == hg002[0]:
            label = "HG002.M:"
        elif length == hg002[1]:
            label = "HG002.P:"
        else:
            label = ""
        label += f"AC={acs[length]}"
        if length == ref_len:
            wrote_ref = True
            label += "_REF"
        ret += f">{label}\n{deduped[length]}\n"
    if not wrote_ref:
        ret += f">REF\n{refseq}\n"
    return ret


def find_name(data, chrom, start, end):
    row = data.loc[chrom, start, end]
    patho = row['patho']
    codis = row['codis']
    name = f'patho_{patho}' if patho != '.' else f'codis_{codis}'

    annos = json.loads(row['annos'])
    motifs = [_['motif'] for _ in annos]
    return name, motifs

if __name__ == '__main__':
    haps = joblib.load("haps.jl")
    # copied from extract haps' temp file
    ref = pysam.FastaFile("refhaps.fa")
    tr_visualizer = TandemRepeatVizWorker()
    catalog = pd.read_csv("/users/u233287/scratch/code/adotto/regions/adotto_TRregions_v1.2.bed", sep='\t')
    sites = (catalog['patho'] != '.') | (catalog['codis'] != '.')
    sites = (sites 
            #& ((catalog['end'] - catalog['start']) <= 1000) 
            & (catalog['patho'] == 'DRD4'))
    data = catalog[sites].set_index(['chr', 'start', 'end'])
    
    #data = pd.read_csv("toextract.bed", sep='\t').set_index(['chr', 'start', 'end'])

    for m_key in haps.keys():
        chrom, start, end = re.split(':|-', m_key)
        start = int(start)
        end = int(end)
        rseq = ref.fetch(m_key)
        m_seqs = dedup(haps[m_key], rseq)
        with open('test.fa', 'w') as fout:
            fout.write(m_seqs)

        try:
            tr_id, motifs = find_name(data, chrom, start, end)
        except Exception:
            continue
        sample_ids, tr_sequences = get_sample_and_sequence_from_fasta("test.fa")
        try:
            tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motifs)
        except Exception:
            print('failed', tr_id)

