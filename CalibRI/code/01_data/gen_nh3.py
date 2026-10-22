# Person-disjoint NHANES arms: participants split into halves A/B by SEQN (same split for all analytes).
# Each half has its own healthy pool (truth = that half's healthy 2.5/97.5 percentiles) and its own excluded pool.
import pandas as pd, numpy as np
rng=np.random.default_rng(20261022)
R='work/nht_raw/'
thy=pd.concat([pd.read_sas(R+f'THYROD_{s}.xpt',format='xport') for s in 'EFG'])
a=pd.read_csv('work/nhtsh/nh_tsh_adults.csv')[['SEQN','healthy','RIDEXPRG']]
thy=thy.merge(a,on='SEQN'); thy=thy[thy.RIDEXPRG!=1]
k=pd.read_csv('data/nhanes/nhanes_hba1c_2011_2023.csv.gz')
k['healthy']=(k.diabetes_dx==0)&(k.hypertension_dx==0)&(k.ckd_kdigo==0)&k.age.between(20,59)
k['path']=k.diabetes_dx==1
allseq=np.union1d(thy.SEQN.unique(),k.SEQN.unique())
perm=rng.permutation(allseq); halfA=set(perm[:len(perm)//2])
def half(s): return np.where(pd.Series(s).isin(halfA),'A','B')
thy['half']=half(thy.SEQN); k['half']=half(k.SEQN)
arms={'TSH':(thy,'LBXTSH1'),'TT4':(thy,'LBXTT4'),'TT3':(thy,'LBXTT3'),'Tg':(thy,'LBXTGN'),'HbA1c':(k,'hba1c_si')}
NS=[300,500,1000]; FR=[0.1,0.3]; REPS=80
meta=[]; vals=[]; i=300000
for an,(df,c) in arms.items():
  for hf in 'AB':
    s=df[df.half==hf]
    h=s[s.healthy][c].dropna().values
    p=(s[s.path] if an=='HbA1c' else s[~s.healthy])[c].dropna().values
    tLL,tUL=np.quantile(h,[.025,.975]); print(an,hf,'healthy',len(h),'path',len(p),'truth',round(tLL,3),round(tUL,3))
    for f in FR:
      for N in NS:
        for r in range(1,REPS+1):
          kk=int(round(N*f)); v=np.concatenate([rng.choice(h,N-kk,replace=False),rng.choice(p,kk,replace=False)]); rng.shuffle(v); i+=1
          meta.append(dict(ds_id=i,source='nhanes',analyte=an,half=hf,n=N,n_nominal=N,frac=f,rep=r,spec_index=-1,trueLL=tLL,trueUL=tUL,decimals=-1))
          vals.append(pd.DataFrame({'ds_id':i,'value':v}))
pd.DataFrame(meta).to_csv('work/ribundle3/meta.csv',index=False)
pd.concat(vals).to_csv('work/ribundle3/values.csv.gz',index=False,compression='gzip')
print('datasets',len(meta))
