"""Bundle 1 = simulated datasets (gen_sim.R) + pooled NHANES datasets (gen_nh.py)."""
import os, pandas as pd
os.makedirs('work/ribundle', exist_ok=True)
pd.concat([pd.read_csv('work/ribuild/sim_meta.csv'), pd.read_csv('work/ribuild/nh_meta.csv')]).to_csv('work/ribundle/meta.csv', index=False)
pd.concat([pd.read_csv('work/ribuild/sim_values.csv'), pd.read_csv('work/ribuild/nh_values.csv')]).to_csv(
    'work/ribundle/values.csv.gz', index=False, compression='gzip')
print('bundle 1 written')
