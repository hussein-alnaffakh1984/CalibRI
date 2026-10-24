# Semi-synthetic real-data datasets from NHANES: healthy pool (known truth) + excluded pool at fixed fractions
import pandas as pd, numpy as np
rng = np.random.default_rng(20261016)
t = pd.read_csv('work/nhtsh/nh_tsh_adults.csv'); t = t[t.RIDEXPRG != 1]
k = pd.read_csv('data/nhanes/nhanes_hba1c_2011_2023.csv.gz')
hk = (k.diabetes_dx==0)&(k.hypertension_dx==0)&(k.ckd_kdigo==0)&k.age.between(20,59)
pools = {
 'TSH':  (t[t.healthy].LBXTSH1.dropna().values, t[~t.healthy].LBXTSH1.dropna().values, 3),
 'FT4':  (t[t.healthy].LBXT4F.dropna().values,  t[~t.healthy].LBXT4F.dropna().values, 2),
 'HbA1c':(k[hk].hba1c_si.dropna().values,        k[k.diabetes_dx==1].hba1c_si.dropna().values, 1),
}
NS=[300,500,1000,2000]; FR=[0.1,0.3]; REPS=200
meta=[]; vals=[]; i=100000
for a,(h,p,dec) in pools.items():
    tLL,tUL = np.quantile(h,[.025,.975])
    print(a,'healthy',len(h),'path',len(p),'truth',round(tLL,3),round(tUL,3))
    for f in FR:
        for N in NS:
            for r in range(1,REPS+1):
                npth=int(round(N*f)); x=np.concatenate([rng.choice(h,N-npth,replace=False),rng.choice(p,npth,replace=False)])
                rng.shuffle(x); i+=1
                meta.append(dict(ds_id=i,source='nhanes',analyte=a,n=N,n_nominal=N,frac=f,rep=r,spec_index=-1,trueLL=tLL,trueUL=tUL,decimals=dec))
                vals.append(pd.DataFrame({'ds_id':i,'value':x}))
pd.DataFrame(meta).to_csv('work/ribuild/nh_meta.csv',index=False); pd.concat(vals).to_csv('work/ribuild/nh_values.csv',index=False)
print('nhanes datasets',len(meta))
