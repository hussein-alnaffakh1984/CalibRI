"""KNHANES subset for refineR (computing time): the first 20 replicates of every cell (240 datasets)."""
import os, pandas as pd
m = pd.read_csv('work/ribundle6/meta.csv'); m = m[m.rep <= 20]
v = pd.read_csv('work/ribundle6/values.csv.gz'); v = v[v.ds_id.isin(m.ds_id)]
os.makedirs('work/ribundle6s', exist_ok=True)
m.to_csv('work/ribundle6s/meta.csv', index=False); v.to_csv('work/ribundle6s/values.csv.gz', index=False, compression='gzip')
print(len(m), 'datasets')
