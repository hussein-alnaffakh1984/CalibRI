# Extra NHANES real arms (2007-2012 thyroid profile): TT4, TT3, Tg. Same healthy/excluded pools as TSH/FT4.
import pandas as pd, numpy as np
rng=np.random.default_rng(20261019)
R='work/nht_raw/'
x=pd.concat([pd.read_sas(R+f'THYROD_{s}.xpt',format='xport') for s in 'EFG'])
a=pd.read_csv('work/nhtsh/nh_tsh_adults.csv')[['SEQN','healthy','RIDEXPRG']]
x=x.merge(a,on='SEQN'); x=x[x.RIDEXPRG!=1]
pools={'TT4':'LBXTT4','TT3':'LBXTT3','Tg':'LBXTGN'}
NS=[300,500,1000,2000]; FR=[0.1,0.3]; REPS=200
meta=[]; vals=[]; i=200000
for an,c in pools.items():
    h=x[x.healthy][c].dropna().values; p=x[~x.healthy][c].dropna().values
    tLL,tUL=np.quantile(h,[.025,.975]); print(an,len(h),len(p),round(tLL,3),round(tUL,3))
    for f in FR:
        for N in NS:
            for r in range(1,REPS+1):
                k=int(round(N*f)); v=np.concatenate([rng.choice(h,N-k,replace=False),rng.choice(p,k,replace=False)]); rng.shuffle(v); i+=1
                meta.append(dict(ds_id=i,source='nhanes',analyte=an,n=N,n_nominal=N,frac=f,rep=r,spec_index=-1,trueLL=tLL,trueUL=tUL,decimals=-1))
                vals.append(pd.DataFrame({'ds_id':i,'value':v}))
pd.DataFrame(meta).to_csv('work/ribundle2/meta.csv',index=False)
pd.concat(vals).to_csv('work/ribundle2/values.csv.gz',index=False,compression='gzip')
print('datasets',len(meta))
