# Table 4, Supplementary Tables S12 (censoring) and S13 (all periods, scheme d vs c) from samawah_conf_d2.py outputs.
import pandas as pd, json
from decimal import Decimal, ROUND_HALF_UP
def r2(x): return str(Decimal(str(float(x))).quantize(Decimal('0.01'), ROUND_HALF_UP))
D = pd.read_csv('/root/samawah/samawah_results_d2.csv'); C = pd.read_csv('/root/samawah/samawah_results_d2c.csv')
E = pd.read_csv('/root/samawah/samawah_excl_results_d2.csv')
for X in (D, C, E): X['period'] = X.period.astype(str)
def get(X, me, lim, sx, pr): return X[(X.method == me) & (X.limit == lim) & (X.sex == sx) & (X.period == pr)].iloc[0]
def cell(r): return f'{r2(r.est)} ({r2(r.cal_lo)}–{r2(r.cal_hi)})'
M = ['reflimR', 'RINet', 'refineR']
t4 = []
for sx, lab in (('F', 'Women'), ('M', 'Men')):
    n = int(get(D, 'reflimR', 'LRL', sx, 'all').n)
    for i, me in enumerate(M):
        t4.append([f'{lab}, n = {n:,}' if i == 0 else '', me, cell(get(D, me, 'LRL', sx, 'all')), cell(get(D, me, 'URL', sx, 'all'))])
s12 = []
for sx, lab in (('F', 'Women'), ('M', 'Men')):
    for i, me in enumerate(M):
        s12.append([lab if i == 0 else '', me, cell(get(D, me, 'LRL', sx, 'all')), cell(get(E, me, 'LRL', sx, 'all')),
                    cell(get(D, me, 'URL', sx, 'all')), cell(get(E, me, 'URL', sx, 'all'))])
names = {'all': '2023–2026', '2324': '2023–24', '2425': '2024–25', '2526': '2025–26'}
s13 = []
for sx, lab in (('F', 'Women'), ('M', 'Men')):
    for pr in ('all', '2324', '2425', '2526'):
        r0 = get(D, 'reflimR', 'LRL', sx, pr)
        for i, me in enumerate(M):
            a, b = get(D, me, 'LRL', sx, pr), get(D, me, 'URL', sx, pr); ca, cb = get(C, me, 'LRL', sx, pr), get(C, me, 'URL', sx, pr)
            s13.append([f'{lab} {names[pr]}' if i == 0 else '', f'{int(r0.n):,} / {int(r0.bin):,}' if i == 0 else '', me,
                        cell(a), f'{r2(ca.cal_lo)}–{r2(ca.cal_hi)}', cell(b), f'{r2(cb.cal_lo)}–{r2(cb.cal_hi)}'])
json.dump({'t4': t4, 's12': s12, 's13': s13}, open('/root/paper13/samawah_tables_d2.json', 'w'), ensure_ascii=False)
full = D[D.period == 'all']; print('full-period URL intervals contain 5.10:', bool(((full[full.limit == 'URL'].cal_lo <= 5.10) & (full[full.limit == 'URL'].cal_hi >= 5.10)).all()))
y = D[(D.period != 'all') & (D.limit == 'URL')]; print('yearly URL intervals excluding 5.10:'); print(y[~((y.cal_lo <= 5.10) & (y.cal_hi >= 5.10))][['method', 'sex', 'period', 'cal_lo', 'cal_hi']])
dc = D.merge(C, on=['method', 'limit', 'sex', 'period'], suffixes=('', '_c')); f = dc[dc.period == 'all']
print('max |d-c| full period bounds:', round(max((f.cal_lo-f.cal_lo_c).abs().max(), (f.cal_hi-f.cal_hi_c).abs().max()), 3))
de = D.merge(E, on=['method', 'limit', 'sex', 'period'], suffixes=('', '_e')); f = de[de.period == 'all']
print('censoring: max |est diff| LRL', round((f[f.limit == 'LRL'].est-f[f.limit == 'LRL'].est_e).abs().max(), 3), 'URL', round((f[f.limit == 'URL'].est-f[f.limit == 'URL'].est_e).abs().max(), 3))
