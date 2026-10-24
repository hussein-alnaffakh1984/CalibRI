# Conformal recalibration of the Al-Samawah TSH limits with the validated scheme (d) composition:
# simulated datasets (bundle 1) + person-disjoint NHANES datasets of the OTHER analytes (total T4, total T3,
# thyroglobulin, HbA1c; bundle 3 for n = 300-1,000, bundle 7 for n = 2,000). Both halves are used, because no
# Al-Samawah patient is an NHANES participant. FT4 (coarse) is never used. Scheme (c) sensitivity: TSH halves added.
# Usage: python samawah_conf_d2.py <scheme d|c> <samawah dir> <out csv>
import sys, pandas as pd, numpy as np
SCHEME, SAM, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
import os; R = os.environ.get('CALIBRI_WORK', 'work/')
def qhat(v, a=0.05):
    v = np.sort(v); k = int(np.ceil((len(v) + 1) * (1 - a))); return v[min(k, len(v)) - 1]
# simulated (bundle 1)
m1 = pd.read_csv(R+'ribundle/meta.csv'); sim = m1[m1.source == 'sim']
fA = pd.concat([pd.read_csv(R+p) for p in ['res_A1/res_fast_shard01.csv', 'res_A2/res_fast_shard02.csv']])
sim = sim.merge(fA, on='ds_id').merge(pd.read_csv(R+'res_C/res_rinet.csv'), on='ds_id')
# person-disjoint NHANES, n = 300-1,000 (bundle 3)
b3 = pd.read_csv(R+'ribundle3/meta.csv').merge(pd.read_csv(R+'res_b3/res_fast_shard01.csv'), on='ds_id').merge(pd.read_csv(R+'res_b3/res_rinet.csv'), on='ds_id')
# person-disjoint NHANES, n = 2,000 (bundle 7): reflimR + RINet (B = 50 from the seeded run) + refineR (separate run)
b7 = pd.read_csv(R+'ribundle7/meta.csv').merge(pd.read_csv(R+'res_b7/reflim_b7.csv'), on='ds_id').merge(pd.read_csv(R+'res_b7/rinet_b7.csv'), on='ds_id')
b7 = b7.rename(columns={'B50_LL_lo': 'rn_CI95_LL_lo', 'B50_LL_hi': 'rn_CI95_LL_hi', 'B50_UL_lo': 'rn_CI95_UL_lo', 'B50_UL_hi': 'rn_CI95_UL_hi'})
rf7 = []
for d in ('res_b7o', 'res_b7t'):
    try: rf7.append(pd.read_csv(R+d+'/res_fast_shard01.csv')[['ds_id','rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi']])
    except FileNotFoundError: pass
if rf7: b7 = b7.merge(pd.concat(rf7), on='ds_id', how='left')
real = pd.concat([b3, b7], ignore_index=True)
keep = ['TT4', 'TT3', 'Tg', 'HbA1c'] + (['TSH'] if SCHEME == 'c' else [])
real = real[real.analyte.isin(keep)]
cal = pd.concat([sim, real], ignore_index=True)
sm = pd.read_csv(SAM+'/meta.csv').merge(pd.read_csv(SAM+'/out/res_fast_shard01.csv'), on='ds_id').merge(pd.read_csv(SAM+'/out/res_rinet.csv'), on='ds_id')
M = {'refineR': ('rf_LL', 'rf_UL', 'rf_UM95_LL_lo', 'rf_UM95_LL_hi', 'rf_UM95_UL_lo', 'rf_UM95_UL_hi'),
     'reflimR': ('rl_LL', 'rl_UL', 'rl_CI95_LL_lo', 'rl_CI95_LL_hi', 'rl_CI95_UL_lo', 'rl_CI95_UL_hi'),
     'RINet': ('rn_LL', 'rn_UL', 'rn_CI95_LL_lo', 'rn_CI95_LL_hi', 'rn_CI95_UL_lo', 'rn_CI95_UL_hi')}
BINS = np.array([300, 500, 1000, 2000]); out = []
for meth, (pL, pU, lLo, lHi, uLo, uHi) in M.items():
    for lim, p, lo, hi, t in (('LRL', pL, lLo, lHi, 'trueLL'), ('URL', pU, uLo, uHi, 'trueUL')):
        c = cal[['source', 'n_nominal', p, lo, hi, t]].dropna(); c = c[(c[p] > 0) & (c[lo] > 0) & (c[hi] > 0)]
        s = np.maximum((np.log(c[hi]) - np.log(c[lo])) / (2 * 1.96), 1e-4); c = c.assign(sc=np.abs(np.log(c[p]) - np.log(c[t])) / s)
        for _, row in sm.iterrows():
            b = BINS[np.argmin(np.abs(np.log(BINS) - np.log(row.n)))]; cb = c[c.n_nominal == b]
            n_sim = int((cb.source == 'sim').sum()); n_real = len(cb) - n_sim
            q = qhat(cb.sc.values) if len(cb) else np.nan
            ss = max((np.log(row[hi]) - np.log(row[lo])) / (2 * 1.96), 1e-4) if row[lo] > 0 else np.nan
            out.append(dict(method=meth, limit=lim, sex=row.sex, period=row.period, n=row.n, bin=b, n_cal_sim=n_sim, n_cal_real=n_real,
                            est=row[p], native_lo=row[lo], native_hi=row[hi], q=q, cal_lo=row[p]*np.exp(-q*ss), cal_hi=row[p]*np.exp(q*ss)))
O = pd.DataFrame(out); O.to_csv(OUT, index=False)
pd.set_option('display.width', 250); print(O.round(3).to_string(index=False))
