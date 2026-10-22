# Supplementary Table S1: simulated benchmark per analyte x n x limit: native coverage, median relative bias, failure rate.
import pandas as pd, numpy as np
m=pd.read_csv('work/ribundle/meta.csv'); m=m[m.source=='sim']
f=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv']]); r=pd.read_csv('work/res_C/res_rinet.csv')
d=m.merge(f,on='ds_id').merge(r,on='ds_id')
M={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','rn_CI95_LL_lo','rn_CI95_LL_hi','rn_CI95_UL_lo','rn_CI95_UL_hi'),
   'refineR':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi')}
rows=[]
for (a,n),g in d.groupby(['analyte','n_nominal']):
  for lim in ('LRL','URL'):
    row=dict(analyte=a,n=n,limit=lim)
    for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
      p,lo,hi,t=(pL,lLo,lHi,'trueLL') if lim=='LRL' else (pU,uLo,uHi,'trueUL')
      ok=g[p].notna()&g[lo].notna()&g[hi].notna()&(g[p]>0)&(g[lo]>0)&(g[hi]>0)
      gg=g[ok]
      row[f'{meth}_cov']=((gg[lo]<=gg[t])&(gg[t]<=gg[hi])).mean()
      row[f'{meth}_bias']=100*np.median(gg[p]/gg[t]-1)
      row[f'{meth}_fail']=100*(1-ok.mean())
    rows.append(row)
S=pd.DataFrame(rows); S.to_csv('work/ana/supp_sim.csv',index=False)
print(S.round(2).to_string(index=False))
