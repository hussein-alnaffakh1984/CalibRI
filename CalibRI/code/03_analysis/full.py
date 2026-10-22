import pandas as pd, numpy as np
m=pd.concat([pd.read_csv('work/ribundle/meta.csv'),pd.read_csv('work/ribundle2/meta.csv')])
f=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv','work/res_b2/res_fast_shard01.csv']])
r=pd.concat([pd.read_csv('work/res_C/res_rinet.csv'),pd.read_csv('work/res_b2/res_rinet.csv')])
d=m.merge(f,on='ds_id',how='left').merge(r,on='ds_id',how='left')
d['arm']=d.source+':'+d.analyte
# reporting resolution tolerance (half step) for coarse real analytes
d['tol']=0.0; d.loc[d.arm=='nhanes:FT4','tol']=0.05; d.loc[d.arm=='nhanes:HbA1c','tol']=0.05; d.loc[d.arm=='nhanes:TT4','tol']=0.05; d.loc[d.arm=='nhanes:TT3','tol']=0.5
print('datasets',len(d), d.groupby('arm').size().to_dict())
M={'refineR':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi'),
   'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','rn_CI95_LL_lo','rn_CI95_LL_hi','rn_CI95_UL_lo','rn_CI95_UL_hi')}
def qhat(s,a=0.05):
    s=np.sort(s); m=len(s); k=int(np.ceil((m+1)*(1-a))); return s[min(k,m)-1]
nat=[]; conf=[]
REAL=['nhanes:TSH','nhanes:HbA1c','nhanes:TT4','nhanes:TT3','nhanes:Tg']
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['arm','source','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['arm','source','n','tol','e','lo','hi','T']
    x['ok']=x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)
    x['cov']=((x.lo<=x['T']+x.tol)&(x['T']-x.tol<=x.hi)).where(x.ok)
    g=x.groupby(['arm','n']).agg(fail=('ok',lambda s:1-s.mean()),cov=('cov','mean'),k=('cov','count')).reset_index(); g['method']=meth; g['limit']=lim; nat.append(g)
    y=x[x.ok].copy(); y['lest']=np.log(y.e); y['s']=np.maximum((np.log(y.hi)-np.log(y.lo))/(2*1.96),1e-4)
    y['nrm']=np.abs(y.lest-np.log(y['T']))/y.s
    pool=y[y.arm!='nhanes:FT4']
    for a in [a for a in pool.arm.unique()]:
      for scheme in ('sim-only','sim+other-real'):
        if scheme=='sim-only' and not a.startswith('nhanes'): continue
        cal=pool[(pool.source=='sim')] if scheme=='sim-only' else pool[pool.arm!=a]
        if scheme=='sim-only': pass
        elif not a.startswith('nhanes'): cal=pool[pool.arm!=a]
        for n in sorted(pool.n.unique()):
          c=cal[cal.n==n]; te=pool[(pool.arm==a)&(pool.n==n)]
          if a.startswith('sim') and scheme=='sim-only': continue
          qq=qhat(c.nrm.values); half=qq*te.s
          cv=((np.exp(te.lest-half)<=te['T']+te.tol)&(te['T']-te.tol<=np.exp(te.lest+half))).mean()
          conf.append(dict(method=meth,limit=lim,arm=a,scheme=scheme,n=n,coverage=cv,hw=np.median(np.exp(half)-1)*100,k=len(te)))
N=pd.concat(nat); C=pd.DataFrame(conf)
N.to_csv('work/ana/final_native.csv',index=False); C.to_csv('work/ana/final_conformal.csv',index=False)
pd.set_option('display.width',250)
print('\n=== NATIVE coverage (new arms) ===')
print(N[N.arm.isin(['nhanes:TT4','nhanes:TT3','nhanes:Tg'])].pivot_table(index=['arm','n'],columns=['method','limit'],values='cov').round(2).to_string())
print('\nfail rates new arms:'); print(N[N.arm.isin(['nhanes:TT4','nhanes:TT3','nhanes:Tg'])&(N.limit=='URL')].groupby(['method','arm']).fail.mean().round(3).unstack(0))
for sc in ('sim-only','sim+other-real'):
  c=C[C.scheme==sc].copy(); c['cell']=c.coverage.round(2).astype(str)+' ('+c.hw.round(0).astype(int).astype(str)+')'
  print(f'\n=== CONFORMAL [{sc}] coverage (median half-width %) ===')
  print(c.pivot_table(index=['arm','n'],columns=['method','limit'],values='cell',aggfunc='first').to_string())
# summary: share of real-arm cells within [0.90,1] 
for sc in ('sim-only','sim+other-real'):
  c=C[(C.scheme==sc)&C.arm.isin(REAL)]; print(sc,'real cells >=0.90:',round((c.coverage>=0.90).mean(),2),'>=0.85:',round((c.coverage>=0.85).mean(),2),'min',round(c.coverage.min(),2))
nr=N[N.arm.isin(REAL)]; print('native real cells >=0.90:',round((nr['cov']>=0.90).mean(),2),'min',round(nr['cov'].min(),2))
