"""NHANES 2007-2012 thyroid analytes: adults (>=20 y) with a TSH result and the healthy-pool flag.

Inputs (CDC public-use files, see data/nhanes/download_nhanes.sh): DEMO, MCQ, RXQ_RX and THYROD,
cycles E (2007-08), F (2009-10), G (2011-12), in work/nht_raw/.
Healthy pool: no thyroid medication (levothyroxine, liothyronine, thyroid, methimazole, propylthiouracil,
liotrix, amiodarone, lithium), no reported thyroid problem (MCQ160M/MCQ170M == 1), TPOAb < 9 IU/mL and
TgAb < 4 IU/mL (missing antibody -> not healthy), not pregnant. 'Don't know' answers and undetermined
pregnancy status are retained. Pregnant participants are dropped later by the dataset generators.
Output: work/nhtsh/nh_tsh_adults.csv
"""
import os, pandas as pd
R = 'work/nht_raw/'
DRUGS = ['LEVOTHYROXINE', 'LIOTHYRONINE', 'THYROID', 'METHIMAZOLE', 'PROPYLTHIOURACIL', 'LIOTRIX', 'AMIODARONE', 'LITHIUM']
out = []
for s in 'EFG':
    th = pd.read_sas(R + f'THYROD_{s}.xpt', format='xport')
    de = pd.read_sas(R + f'DEMO_{s}.xpt', format='xport')[['SEQN', 'RIAGENDR', 'RIDAGEYR', 'RIDEXPRG']]
    mc = pd.read_sas(R + f'MCQ_{s}.xpt', format='xport')[['SEQN', 'MCQ160M', 'MCQ170M']]
    rx = pd.read_sas(R + f'RXQ_RX_{s}.xpt', format='xport', encoding='latin1')
    d = rx.RXDDRUG.fillna('').str.upper()
    med = set(rx.SEQN[d.apply(lambda x: any(k in x.split('; ') or k == x or x.startswith(k) for k in DRUGS))])
    m = th[['SEQN', 'LBXTSH1', 'LBXT4F', 'LBXTPO', 'LBXATG']].merge(de, on='SEQN').merge(mc, on='SEQN', how='left')
    m = m[(m.RIDAGEYR >= 20) & m.LBXTSH1.notna()].copy()
    m['cycle'] = s; m['thymed'] = m.SEQN.isin(med)
    out.append(m)
m = pd.concat(out, ignore_index=True)
m['healthy'] = (~m.thymed) & (m.MCQ160M != 1) & (m.MCQ170M != 1) & (m.LBXTPO < 9) & (m.LBXATG < 4) & (m.RIDEXPRG != 1)
os.makedirs('work/nhtsh', exist_ok=True)
m.to_csv('work/nhtsh/nh_tsh_adults.csv', index=False)
print('adults', len(m), '| healthy', int(m.healthy.sum()), '| thyroid medication', int(m.thymed.sum()))
