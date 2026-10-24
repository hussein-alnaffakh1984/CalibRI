# Repeated person-disjoint splits: same design as bundle 3 (gen_nh3.py), with 4 additional random splits.
# Split k uses its own seed; each split: 5 analytes x 2 halves x 2 fractions x 3 n x 80 reps = 4,800 datasets.
import pandas as pd, numpy as np
R='work/nht_raw/'
thy=pd.concat([pd.read_sas(R+f'THYROD_{s}.xpt',format='xport') for s in 'EFG'])
a=pd.read_csv('work/nhtsh/nh_tsh_adults.csv')[['SEQN','healthy','RIDEXPRG']]
thy=thy.merge(a,on='SEQN'); thy=thy[thy.RIDEXPRG!=1]
k=pd.read_csv('data/nhanes/nhanes_hba1c_2011_2023.csv.gz')
k['healthy']=(k.diabetes_dx==0)&(k.hypertension_dx==0)&(k.ckd_kdigo==0)&k.age.between(20,59)
k['path']=k.diabetes_dx==1
allseq=np.union1d(thy.SEQN.unique(),k.SEQN.unique())
arms={'TSH':(thy,'LBXTSH1'),'TT4':(thy,'LBXTT4'),'TT3':(thy,'LBXTT3'),'Tg':(thy,'LBXTGN'),'HbA1c':(k,'hba1c_si')}
NS=[300,500,1000]; FR=[0.1,0.3]; REPS=80
meta=[]; vals=[]; i=400000; truth=[]
for split in (2,3,4,5):                       # split 1 = bundle 3 (seed 20261022)
  rng=np.random.default_rng(20261022+split)
  perm=rng.permutation(allseq); halfA=set(perm[:len(perm)//2])
  for an,(df,c) in arms.items():
    hcol=np.where(df.SEQN.isin(halfA),'A','B')
    for hf in 'AB':
      s=df[hcol==hf]
      h=s[s.healthy][c].dropna().values
      p=(s[s.path] if an=='HbA1c' else s[~s.healthy])[c].dropna().values
      tLL,tUL=np.quantile(h,[.025,.975]); truth.append(dict(split=split,analyte=an,half=hf,n_healthy=len(h),n_excl=len(p),trueLL=tLL,trueUL=tUL))
      for f in FR:
        for N in NS:
          for r in range(1,REPS+1):
            kk=int(round(N*f)); v=np.concatenate([rng.choice(h,N-kk,replace=False),rng.choice(p,kk,replace=False)]); rng.shuffle(v); i+=1
            meta.append(dict(ds_id=i,split=split,source='nhanes',analyte=an,half=hf,n=N,n_nominal=N,frac=f,rep=r,trueLL=tLL,trueUL=tUL))
            vals.append(pd.DataFrame({'ds_id':i,'value':v}))
pd.DataFrame(truth).to_csv('work/ribundle4/truth.csv',index=False)
pd.DataFrame(meta).to_csv('work/ribundle4/meta.csv',index=False)
pd.concat(vals).to_csv('work/ribundle4/values.csv.gz',index=False,compression='gzip')
print('datasets',len(meta))
