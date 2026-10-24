# RINet (LeBien, Velev & Roche-Lima, J Biomed Inform 2026) 1D inference with B bootstrap resamples.
# Re-implementation of the rinet 0.1.1 R package interface (bundled model + scaler), vectorised histogram features.
# Reports percentile CIs from the first 50 resamples and from all B resamples. Each dataset has its own seed (= ds_id),
# so results do not depend on run order or restarts.
# Usage: python rinet_boot.py <data_dir> <model_dir> <out_csv> [B=1000] [filter: all|sim_small]
import os, sys, pickle, time, numpy as np, pandas as pd
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
import keras
from scipy.stats import norm
DATA, MODELS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
B = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
FILT = sys.argv[5] if len(sys.argv) > 5 else 'all'
model = keras.models.load_model(os.path.join(MODELS, "rinet_1d.keras"))
scaler = pickle.load(open(os.path.join(MODELS, "scaler_1d.pkl"), "rb"))
BR = np.linspace(-4, 4, 101); Z = norm.ppf([0.025, 0.975])

def feats(S):                                   # identical to np.histogram(density=True) + min-max scaling
    Zs = (S - S.mean(1, keepdims=True)) / S.std(1, keepdims=True)
    idx = np.searchsorted(BR, Zs, side='right') - 1; idx[Zs == BR[-1]] = 99
    ok = (Zs >= BR[0]) & (Zs <= BR[-1])
    H = np.zeros((S.shape[0], 100)); rows = np.repeat(np.arange(S.shape[0]), S.shape[1])
    np.add.at(H, (rows[ok.ravel()], idx.ravel()[ok.ravel()]), 1)
    mn = H.min(1, keepdims=True); mx = H.max(1, keepdims=True); return (H - mn) / (mx - mn)

meta = pd.read_csv(os.path.join(DATA, "meta.csv"))
if FILT == 'sim_small': meta = meta[(meta.source == 'sim') & meta.n_nominal.isin([300, 500, 1000])]
vals = pd.read_csv(os.path.join(DATA, "values.csv.gz"))
vals = vals[vals.ds_id.isin(set(meta.ds_id))]
groups = {k: v.values for k, v in vals.groupby("ds_id")["value"]}
done = set(pd.read_csv(OUT).ds_id) if os.path.exists(OUT) else set()
ids = [i for i in meta.ds_id if i not in done]
print(f"RINet B={B}: {len(ids)} datasets to run ({len(done)} done)", flush=True)
t0 = time.time(); rows = []
for k, i in enumerate(ids, 1):
    x = groups[i].astype(float); nonpos = int((x <= 0).sum())
    if nonpos: x[x <= 0] = x[x > 0].min() / 2
    lx = np.log(x); rng = np.random.default_rng(int(i))
    S = np.vstack([lx[None], lx[rng.integers(0, len(lx), (B, len(lx)))]])
    pred = scaler.inverse_transform(model.predict(feats(S)[..., None], verbose=0, batch_size=B + 1))
    m = S.mean(1); s = S.std(1)
    mu = pred[:, 0] * s + m; sd = pred[:, 1] * s
    LL = np.exp(mu + Z[0] * sd); UL = np.exp(mu + Z[1] * sd)
    r = dict(ds_id=i, rn_LL=LL[0], rn_UL=UL[0], rn_nonpos=nonpos)
    for tag, sl in (('B50', slice(1, 51)), (f'B{B}', slice(1, B + 1))):
        r[f'{tag}_LL_lo'], r[f'{tag}_LL_hi'] = np.quantile(LL[sl], [0.025, 0.975])
        r[f'{tag}_UL_lo'], r[f'{tag}_UL_hi'] = np.quantile(UL[sl], [0.025, 0.975])
    rows.append(r)
    if k % 200 == 0 or k == len(ids):
        pd.DataFrame(rows).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False); rows = []
        print(f"{k}/{len(ids)} ({(time.time()-t0)/60:.1f} min)", flush=True)
print("Finished.")
