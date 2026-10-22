# External validation benchmark from KNHANES 2013-2015 (Korea): TSH and FT4, healthy pool + excluded pool.
# Healthy: age >=20, not pregnant, no physician-diagnosed thyroid disease (DE2_dg) or thyroid cancer (DC7_dg),
# no thyroid treatment (TH_tx), TPOAb <=34 IU/mL (positivity >34.0 IU/mL, Kim et al. Endocrinol Metab 2017).
# Excluded pool: all other adults (non-pregnant) with a result. 2 analytes x 2 fractions x 3 n x 200 reps = 2,400 datasets.
import pyreadstat, pandas as pd, numpy as np
cols=['id','year','sex','age','HE_prg','DE2_dg','DC7_dg','TH_tx','HE_TSH','HE_fT4','HE_TPOAb']
D=pd.concat([pyreadstat.read_sas7bdat(f'work/knhanes/raw/hn{y}_all.sas7bdat',usecols=cols,encoding='cp949')[0] for y in (13,14,15)],ignore_index=True)
D=D[D.HE_TSH.notna()&(D.age>=20)&(D.HE_prg!=1)]
D['healthy']=(D.DE2_dg!=1)&(D.DC7_dg!=1)&(D.TH_tx!=1)&(D.HE_TPOAb<=34)
rng=np.random.default_rng(20261022)
meta=[]; vals=[]; i=600000; summ=[]
for an,c in (('TSH','HE_TSH'),('FT4','HE_fT4')):
    h=D[D.healthy][c].dropna().values; p=D[~D.healthy][c].dropna().values
    tLL,tUL=np.quantile(h,[.025,.975]); summ.append(dict(analyte=an,n_healthy=len(h),n_excl=len(p),trueLL=tLL,trueUL=tUL,
        excl_outside=np.mean((p<tLL)|(p>tUL))))
    for f in (0.1,0.3):
        for N in (300,500,1000):
            for r in range(1,201):
                k=int(round(N*f)); v=np.concatenate([rng.choice(h,N-k,replace=False),rng.choice(p,k,replace=False)]); rng.shuffle(v); i+=1
                meta.append(dict(ds_id=i,source='knhanes',analyte=an,n=N,n_nominal=N,frac=f,rep=r,trueLL=tLL,trueUL=tUL))
                vals.append(pd.DataFrame({'ds_id':i,'value':v}))
pd.DataFrame(summ).to_csv('work/ribundle6/truth.csv',index=False); print(pd.DataFrame(summ).round(3).to_string())
pd.DataFrame(meta).to_csv('work/ribundle6/meta.csv',index=False)
pd.concat(vals).to_csv('work/ribundle6/values.csv.gz',index=False,compression='gzip'); print('datasets',len(meta))
