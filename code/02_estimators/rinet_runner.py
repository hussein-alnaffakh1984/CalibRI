# RINet (LeBien, Velev & Roche-Lima, J Biomed Inform 2026) 1D inference, re-implemented in Python
# from the official R package rinet 0.1.1 (CRAN, MIT licence), using its bundled model + scaler files.
# Usage: python rinet_runner.py <data_dir> <model_dir> <out_dir> [n_boot=50] [max_ds]
import os, sys, pickle, time, numpy as np, pandas as pd
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
import keras
from scipy.stats import norm

DATA, MODELS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NBOOT = int(sys.argv[4]) if len(sys.argv) > 4 else 50
MAXDS = int(sys.argv[5]) if len(sys.argv) > 5 else None
os.makedirs(OUT, exist_ok=True)
model = keras.models.load_model(os.path.join(MODELS, "rinet_1d.keras"))
scaler = pickle.load(open(os.path.join(MODELS, "scaler_1d.pkl"), "rb"))
BREAKS = np.linspace(-4, 4, 101)
Z = norm.ppf([0.025, 0.975])

def features(logx):
    m = logx.mean(); s = logx.std()                      # population SD, as in the R package
    h, _ = np.histogram((logx - m) / s, bins=BREAKS, density=True)
    h = (h - h.min()) / (h.max() - h.min())
    return h, m, s

meta = pd.read_csv(os.path.join(DATA, "meta.csv"))
vals = pd.read_csv(os.path.join(DATA, "values.csv.gz"))
groups = {k: v.values for k, v in vals.groupby("ds_id")["value"]}
ids = meta.ds_id.tolist()[:MAXDS] if MAXDS else meta.ds_id.tolist()
resfile = os.path.join(OUT, "res_rinet.csv")
done = set(pd.read_csv(resfile).ds_id) if os.path.exists(resfile) else set()
ids = [i for i in ids if i not in done]
print(f"RINet: {len(ids)} datasets to run ({len(done)} already done), {NBOOT} bootstraps each", flush=True)
rng = np.random.default_rng(20261016)
t0 = time.time(); rows = []
for k, i in enumerate(ids, 1):
    x = groups[i].astype(float)
    n_nonpos = int((x <= 0).sum())
    if n_nonpos:                                         # log scale needs x > 0: replace by half the smallest positive value
        x[x <= 0] = x[x > 0].min() / 2
    lx = np.log(x)
    samples = [lx] + [rng.choice(lx, len(lx), replace=True) for _ in range(NBOOT)]
    F, M, S = zip(*[features(s) for s in samples])
    pred = scaler.inverse_transform(model.predict(np.stack(F)[..., None], verbose=0))
    mu = pred[:, 0] * np.array(S) + np.array(M); sd = pred[:, 1] * np.array(S)
    LL = np.exp(mu + Z[0] * sd); UL = np.exp(mu + Z[1] * sd)
    r = dict(ds_id=i, rn_LL=LL[0], rn_UL=UL[0], rn_frac=pred[0, 2] if pred.shape[1] > 2 else np.nan, rn_nonpos=n_nonpos)
    if NBOOT:
        for lev, a in ((95, 0.025), (90, 0.05)):
            r[f"rn_CI{lev}_LL_lo"], r[f"rn_CI{lev}_LL_hi"] = np.quantile(LL[1:], [a, 1 - a])
            r[f"rn_CI{lev}_UL_lo"], r[f"rn_CI{lev}_UL_hi"] = np.quantile(UL[1:], [a, 1 - a])
        r["rn_boot_sdlogUL"] = np.std(np.log(UL[1:])); r["rn_boot_sdlogLL"] = np.std(np.log(LL[1:]))
    rows.append(r)
    if k % 200 == 0 or k == len(ids):
        pd.DataFrame(rows).to_csv(resfile, mode="a", header=not os.path.exists(resfile), index=False); rows = []
        print(f"{k}/{len(ids)} done ({(time.time()-t0)/60:.1f} min)", flush=True)
print("Finished RINet.")
