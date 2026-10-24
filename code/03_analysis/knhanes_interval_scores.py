# KNHANES external validation: coverage and median relative interval score (Supplementary Table S9, KNHANES column)
# for native intervals, scheme (d) calibration (simulated + person-disjoint NHANES other analytes) and real-only calibration.
# RINet uses the B = 50 intervals of the per-dataset-seeded run.
import pandas as pd, numpy as np
import os; W=os.environ.get('CALIBRI_WORK','work/')
def qhat(v,a=0.05):
    v=np.sort(np.asarray(v)); m=len(v); k=int(np.ceil((m+1)*(1-a))); return v[min(k,m)-1]
RL=['ds_id','rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi']
m1=pd.read_csv(W+'ribundle/meta.csv'); m1=m1[(m1.source=='sim')&m1.n_nominal.isin([300,500,1000])]
f1=pd.concat([pd.read_csv(W+'res_A1/res_fast_shard01.csv'),pd.read_csv(W+'res_A2/res_fast_shard02.csv')])
sim=m1.merge(f1[RL],on='ds_id').merge(pd.read_csv(W+'res_b4/rinet_sim.csv'),on='ds_id'); sim['src']='sim'
m3=pd.read_csv(W+'ribundle3/meta.csv'); nh=m3.merge(pd.read_csv(W+'res_b3/res_fast_shard01.csv')[RL],on='ds_id').merge(pd.read_csv(W+'res_b4/rinet_b3.csv'),on='ds_id'); nh['src']='nh'
kn=pd.read_csv(W+'ribundle6/meta.csv').merge(pd.read_csv(W+'res_b6/reflim_b6.csv'),on='ds_id').merge(pd.read_csv(W+'res_b6/rinet_b6.csv'),on='ds_id'); kn['src']='kn'
d=pd.concat([sim,nh,kn],ignore_index=True)
V={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),'RINet':('rn_LL','rn_UL','B50_LL_lo','B50_LL_hi','B50_UL_lo','B50_UL_hi')}
rows=[]
for meth,(pL,pU,a1,a2,b1,b2) in V.items():
  for lim,p,lo,hi,t in (('LRL',pL,a1,a2,'trueLL'),('URL',pU,b1,b2,'trueUL')):
    x=d[['src','analyte','n_nominal',p,lo,hi,t]].copy(); x.columns=['src','an','n','e','lo','hi','T']
    x=x.dropna(); x=x[(x.e>0)&(x.lo>0)&(x.hi>0)].copy(); x['s']=np.maximum(np.log(x.hi/x.lo)/3.92,1e-4); x['sc']=np.abs(np.log(x.e/x['T']))/x.s
    for a in ('TSH','FT4'):
      for n in (300,500,1000):
        te=x[(x.src=='kn')&(x.an==a)&(x.n==n)]; S=x[(x.src=='sim')&(x.n==n)]; N=x[(x.src=='nh')&(x.n==n)&(x.an!=a)]
        for sc,c in (('native',None),('d',pd.concat([S,N])),('realonly',N)):
          if c is None: L,U=te.lo,te.hi
          else: q=qhat(c.sc); L,U=te.e*np.exp(-q*te.s),te.e*np.exp(q*te.s)
          T=te['T']; IS=((U-L)+40*np.maximum(L-T,0)+40*np.maximum(T-U,0))/T
          rows.append(dict(meth=meth,lim=lim,a=a,n=n,sc=sc,cov=((L<=T)&(T<=U)).mean(),IS=np.median(IS)))
R=pd.DataFrame(rows)
print(R.groupby('sc').apply(lambda g: pd.Series({'p90':100*(g['cov']>=.9).mean(),'med':g['cov'].median(),'min':g['cov'].min(),'IS':g.IS.median()}),include_groups=False).round(4))
