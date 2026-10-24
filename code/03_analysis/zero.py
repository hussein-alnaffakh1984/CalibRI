# Bias and native coverage at 0% contamination vs mixed datasets (10/30% pooled), split-1 halves, reflimR and RINet.
import pandas as pd, numpy as np
TOL={'TT4':0.05,'TT3':0.5,'HbA1c':0.05}
m5=pd.read_csv('work/ribundle5/meta.csv'); z=m5.merge(pd.read_csv('work/res_b5/reflim_b5.csv'),on='ds_id').merge(pd.read_csv('work/res_b5/rinet_b5.csv'),on='ds_id')
m3=pd.read_csv('work/ribundle3/meta.csv'); f3=pd.read_csv('work/res_b3/res_fast_shard01.csv')[['ds_id','rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi']]
mx=m3.merge(f3,on='ds_id').merge(pd.read_csv('work/res_b4/rinet_b3.csv'),on='ds_id')
V={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','B50_LL_lo','B50_LL_hi','B50_UL_lo','B50_UL_hi')}
rows=[]
for arm,d in (('0%',z),('10–30%',mx)):
  d=d.copy(); d['tol']=d.analyte.map(TOL).fillna(0)
  for meth,(pL,pU,lLo,lHi,uLo,uHi) in V.items():
    for lim,p,lo,hi,t in (('LRL',pL,lLo,lHi,'trueLL'),('URL',pU,uLo,uHi,'trueUL')):
      x=d[d[p].notna()&(d[p]>0)&d[lo].notna()&d[hi].notna()]
      for (a,n,h),g in x.groupby(['analyte','n_nominal','half']):
        rows.append(dict(arm=arm,method=meth,limit=lim,analyte=a,n=n,half=h,bias=100*np.median(g[p]/g[t]-1),
                         cov=((g[lo]<=g[t]+g.tol)&(g[t]-g.tol<=g[hi])).mean()))
R=pd.DataFrame(rows).groupby(['arm','method','analyte','limit','n'])[['bias','cov']].mean().reset_index()
R.to_csv('work/ana/zero_results.csv',index=False)
P=R.pivot_table(index=['method','analyte','limit'],columns=['arm','n'],values='bias').round(1)
print(P.to_string())
C=R.groupby(['arm','method']).apply(lambda g:pd.Series({'p90':100*(g['cov']>=0.9).mean(),'med':g['cov'].median()}),include_groups=False)
print(C.round(3))
