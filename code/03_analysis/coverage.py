import pandas as pd, numpy as np
m=pd.read_csv('work/ribundle/meta.csv')
f=pd.concat([pd.read_csv('work/res_A1/res_fast_shard01.csv'),pd.read_csv('work/res_A2/res_fast_shard02.csv')])
r=pd.read_csv('work/res_C/res_rinet.csv')
d=m.merge(f,on='ds_id',how='left').merge(r,on='ds_id',how='left')
print('fast rows',len(f),'unique',f.ds_id.nunique(),'| merged',len(d))
d['tol']=np.where((d.source=='nhanes')&d.analyte.isin(['FT4','HbA1c']),0.05,0.0)
M={'refineR':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi'),
   'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','rn_CI95_LL_lo','rn_CI95_LL_hi','rn_CI95_UL_lo','rn_CI95_UL_hi')}
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    ok=d[[p,lo,hi]].notna().all(axis=1)
    cov=((d[lo]<=d[t]+d.tol)&(d[t]-d.tol<=d[hi]))
    tmp=pd.DataFrame({'source':d.source,'analyte':d.analyte,'n':d.n_nominal,'ok':ok,'cov':cov.where(ok),'bias':((d[p]/d[t]-1)*100).where(ok)})
    g=tmp.groupby(['source','analyte','n']).agg(fail=('ok',lambda x:1-x.mean()),cov=('cov','mean'),bias=('bias','median')).reset_index()
    g['method']=meth; g['limit']=lim; rows.append(g)
G=pd.concat(rows); G.to_csv('work/ana/coverage_native.csv',index=False)
pd.set_option('display.width',250)
for lim in ('URL','LRL'):
  t=G[G.limit==lim].pivot_table(index=['source','analyte','n'],columns='method',values='cov').round(2)
  print(f'\n=== native 95% interval coverage, {lim} ==='); print(t.to_string())
print('\nfailure rates:'); print(G[G.limit=='URL'].groupby(['method','source','analyte']).fail.mean().round(3).unstack(0).to_string())
d.to_pickle('work/ana/all.pkl')
