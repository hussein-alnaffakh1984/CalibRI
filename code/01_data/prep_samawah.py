"""Build the Samawah TSH datasets (code only; the laboratory exports are not public).

Input: six Excel exports from the Mindray CL-900i, named {F|M}{2324|2425|2526}.xlsx
       (sex, export period 15 Sep–14 Sep). Only rows with Chemistry == 'TSH' are used.
Cleaning: consecutive exports overlap on 15 September, so a result appearing in an earlier
export (same sex, run date-time and result) is dropped from the later one (8 records);
records without a result (7) are removed. Censored results '<0.005' and '>100' are set to
0.0025 and 100 mIU/L (main analysis) or removed (sensitivity analysis, --exclude-censored).
Output: <out>/meta.csv (ds_id, sex, period, n) and <out>/values.csv.gz (ds_id, value).
Usage:  python prep_samawah.py <excel_dir> <out_dir> [--exclude-censored]
"""
import sys, os, pandas as pd

def load(d):
    rows = []
    for sex in 'FM':
        for per in ('2324', '2425', '2526'):
            x = pd.read_excel(os.path.join(d, f'{sex}{per}.xlsx'), dtype=str)
            x = x[x['Chemistry'].str.strip() == 'TSH']
            rows.append(pd.DataFrame({'sex': sex, 'file': int(per),
                                      'raw': x['Result'].str.strip().values, 'run': x['Run Date'].values}))
    t = pd.concat(rows, ignore_index=True)
    dup = t.duplicated(subset=['sex', 'run', 'raw'], keep='first') & t.run.notna()
    t = t[~dup].copy()
    t['censored'] = t.raw.str.startswith('<') | t.raw.str.startswith('>')
    t['v'] = pd.to_numeric(t.raw.replace({'<0.005': '0.0025', '>100': '100'}), errors='coerce')
    print(f'duplicates removed: {dup.sum()}, empty results removed: {t.v.isna().sum()}')
    return t[t.v.notna()]

if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]; excl = '--exclude-censored' in sys.argv
    t = load(src)
    if excl: t = t[~t.censored]
    os.makedirs(out, exist_ok=True)
    meta, vals, i = [], [], 0
    for sex, per in [('F', 'all'), ('M', 'all')] + [(x, p) for x in 'FM' for p in ('2324', '2425', '2526')]:
            s = t[t.sex == sex] if per == 'all' else t[(t.sex == sex) & (t.file == int(per))]
            i += 1; meta.append(dict(ds_id=i, sex=sex, period=per, n=len(s)))
            vals.append(pd.DataFrame({'ds_id': i, 'value': s.v.values}))
    pd.DataFrame(meta).to_csv(os.path.join(out, 'meta.csv'), index=False)
    pd.concat(vals).to_csv(os.path.join(out, 'values.csv.gz'), index=False, compression='gzip')
    print(pd.DataFrame(meta).to_string(index=False))
