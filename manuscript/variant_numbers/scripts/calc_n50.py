import sys
import truvari

def fasta_reader(fname):
    fh = truvari.opt_gz_open(fname)
    all_lens = []
    cur_len = 0
    next(fh) # Skip the first header
    for line in fh:
        if line.startswith(">"):
            all_lens.append(cur_len)
            cur_len = 0
            continue
        cur_len += len(line.strip())
    all_lens.append(cur_len)
    return all_lens

def calc_n50(lens):
    lens.sort()
    tot_len = sum(lens)
    thresh_len = tot_len / 2
    obs_len = 0
    cur_pos = 0
    while obs_len < thresh_len:
        obs_len += lens[cur_pos]
        cur_pos += 1
    return lens[cur_pos]

if __name__ == '__main__':
    fn = sys.argv[1]
    all_lens = fasta_reader(fn)
    print(len(all_lens), calc_n50(all_lens))
