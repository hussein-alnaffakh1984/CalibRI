# refineR bootstrap confidence interval (50 resamples): coverage of the true TSH limits at n = 500
# (100 simulated, 40 pooled NHANES datasets), with Wilson 95% CIs (main text, Native coverage).
import os, pandas as pd, numpy as np
from scipy.stats import norm
W = os.environ.get('CALIBRI_WORK', 'work/')
b = pd.read_csv(W+'res_B/B_merged.csv')[['ds_id', 'rfb_CI95_LL_lo', 'rfb_CI95_LL_hi', 'rfb_CI95_UL_lo', 'rfb_CI95_UL_hi']].merge(pd.read_csv(W+'ribundle/meta.csv'), on='ds_id')
def wilson(k, n, z=norm.ppf(0.975)):
    p = k/n; c = (p + z*z/(2*n))/(1 + z*z/n); h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1 + z*z/n); return c-h, c+h
for src, g in b.groupby('source'):
    for lim, lo, hi, t in (('URL', 'rfb_CI95_UL_lo', 'rfb_CI95_UL_hi', 'trueUL'), ('LRL', 'rfb_CI95_LL_lo', 'rfb_CI95_LL_hi', 'trueLL')):
        k = int(((g[lo] <= g[t]) & (g[t] <= g[hi])).sum()); n = len(g); a, z = wilson(k, n)
        print(f'{src} {lim}: {k}/{n} = {k/n:.2f} (95% CI {a:.2f}–{z:.2f})')
