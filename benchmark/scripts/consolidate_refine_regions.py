import sys
import pandas as pd

ei = pd.read_csv(sys.argv[1], sep='\t').set_index(["chrom", "start", "end"])
li = pd.read_csv(sys.argv[2], sep='\t').set_index(["chrom", "start", "end"])
th = pd.read_csv(sys.argv[3], sep='\t').set_index(["chrom", "start", "end"])

th['ei'] = ei['state']
th['li'] = li['state']
th[['ei', 'li', 'state']].reset_index().to_csv("/dev/stdout", sep='\t', header=False, index=False)
