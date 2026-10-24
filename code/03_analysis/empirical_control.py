# Empirical-percentile control for the healthy-only (0% contamination) arm (Supplementary Table S7, last column):
# median relative bias (%) of the 2.5th/97.5th sample percentiles of each healthy-only dataset, split 1, by analyte, limit and n.
import os, pandas as pd, numpy as np
W = os.environ.get('CALIBRI_WORK', 'work/')
m = pd.read_csv(W+'ribundle5/meta.csv'); v = pd.read_csv(W+'ribundle5/values.csv.gz')
q = v.groupby('ds_id').value.quantile([0.025, 0.975]).unstack(); q.columns = ['eLL', 'eUL']
d = m.merge(q, left_on='ds_id', right_index=True)
d['bLL'] = 100*(d.eLL/d.trueLL-1); d['bUL'] = 100*(d.eUL/d.trueUL-1)
out = d.groupby(['analyte', 'n_nominal', 'half'])[['bLL', 'bUL']].median().groupby(level=[0, 1]).mean()
print(out.round(1).to_string()); print('max |bias| %:', round(out.abs().max().max(), 2))
