# 0% contamination arm: healthy-pool-only datasets on the split-1 halves (same split as gen_nh3.py, seed 20261022).
# Separates estimator/shape bias from contamination bias. 5 analytes x 2 halves x 3 n x 80 reps = 2,400 datasets.
import pandas as pd, numpy as np
rng=np.random.default_rng(20261022)
R='work/nht_raw/'
thy=pd.concat([pd.read_sas(R+f'THYROD_{s}.xpt',format='xport') for s in 'EFG'])
a=pd.read_csv('work/nhtsh/nh_tsh_adults.csv')[['SEQN','healthy','RIDEXPRG']]
thy=thy.merge(a,on='SEQN'); thy=thy[thy.RIDEXPRG!=1]
k=pd.read_csv('data/nhanes/nhanes_hba1c_2011_2023.csv.gz')
k['healthy']=(k.diabetes_dx==0)&(k.hypertension_dx==0)&(k.ckd_kdigo==0)&k.age.between(20,59)
allseq=np.union1d(thy.SEQN.unique(),k.SEQN.unique())
perm=rng.permutation(allseq); halfA=set(perm[:len(perm)//2])     # identical split to bundle 3
rng2=np.random.default_rng(20261023)
arms={'TSH':(thy,'LBXTSH1'),'TT4':(thy,'LBXTT4'),'TT3':(thy,'LBXTT3'),'Tg':(thy,'LBXTGN'),'HbA1c':(k,'hba1c_si')}
m3=pd.read_csv('work/ribundle3/meta.csv')
meta=[]; vals=[]; i=500000
for an,(df,c) in arms.items():
  for hf in 'AB':
    s=df[np.where(df.SEQN.isin(halfA),'A','B')==hf]; h=s[s.healthy][c].dropna().values
    tLL,tUL=np.quantile(h,[.025,.975])
    ref=m3[(m3.analyte==an)&(m3.half==hf)].iloc[0]; assert abs(ref.trueUL-tUL)<1e-9 and abs(ref.trueLL-tLL)<1e-9
    for N in (300,500,1000):
      for r in range(1,81):
        v=rng2.choice(h,N,replace=False); i+=1
        meta.append(dict(ds_id=i,source='nhanes',analyte=an,half=hf,n=N,n_nominal=N,frac=0.0,rep=r,trueLL=tLL,trueUL=tUL))
        vals.append(pd.DataFrame({'ds_id':i,'value':v}))
pd.DataFrame(meta).to_csv('work/ribundle5/meta.csv',index=False)
pd.concat(vals).to_csv('work/ribundle5/values.csv.gz',index=False,compression='gzip')
print('datasets',len(meta),'(truth matches bundle 3)')
