# External validation on KNHANES 2013-2015 (Korea): calibration uses ONLY simulated + NHANES data.
import pandas as pd, numpy as np, os
def qhat(v,a=0.05):
    v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
RL=['ds_id','rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi']
RF=['ds_id','rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi']
m1=pd.read_csv('work/ribundle/meta.csv'); m1=m1[(m1.source=='sim')&m1.n_nominal.isin([300,500,1000])]
f1=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv']])
sim=m1.merge(f1[RL+RF[1:]],on='ds_id').merge(pd.read_csv('work/res_b4/rinet_sim.csv'),on='ds_id'); sim['src']='sim'
m3=pd.read_csv('work/ribundle3/meta.csv'); f3=pd.read_csv('work/res_b3/res_fast_shard01.csv')
nh=m3.merge(f3[RL+RF[1:]],on='ds_id').merge(pd.read_csv('work/res_b4/rinet_b3.csv'),on='ds_id'); nh['src']='nhanes'
m6=pd.read_csv('work/ribundle6/meta.csv'); kn=m6.merge(pd.read_csv('work/res_b6/reflim_b6.csv'),on='ds_id').merge(pd.read_csv('work/res_b6/rinet_b6.csv'),on='ds_id')
rfk=None
if os.path.exists('work/res_b6/refine/res_fast_shard01.csv'):
    rfk=pd.read_csv('work/res_b6/refine/res_fast_shard01.csv')[RF]; kn=kn.merge(rfk,on='ds_id',how='left')
kn['src']='knhanes'
d=pd.concat([sim,nh,kn],ignore_index=True)
V={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet_B50':('rn_LL','rn_UL','B50_LL_lo','B50_LL_hi','B50_UL_lo','B50_UL_hi'),
   'RINet_B1000':('rn_LL','rn_UL','B1000_LL_lo','B1000_LL_hi','B1000_UL_lo','B1000_UL_hi'),
   'refineR_UM':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi')}
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in V.items():
  if pL not in d: continue
  for lim,p,lo,hi,t in (('LRL',pL,lLo,lHi,'trueLL'),('URL',pU,uLo,uHi,'trueUL')):
    x=d[['src','analyte','n_nominal',p,lo,hi,t]].copy(); x.columns=['src','analyte','n','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(x.lest-np.log(x['T']))/x.s
    for a in ('TSH','FT4'):
      for n in (300,500,1000):
        te=x[(x.src=='knhanes')&(x.analyte==a)&(x.n==n)]
        if len(te)==0: continue
        S=x[(x.src=='sim')&(x.n==n)]; N=x[(x.src=='nhanes')&(x.n==n)]
        for sc,c in (('native',None),('sim-only',S),('sim+NHANES(all)',pd.concat([S,N])),('sim+NHANES(other)',pd.concat([S,N[N.analyte!=a]]))):
          if c is None: L,U=te.lo,te.hi
          else: q=qhat(c.sc.values); L,U=np.exp(te.lest-q*te.s),np.exp(te.lest+q*te.s)
          T=te['T']; cov=((L<=T)&(T<=U)).mean()
          rows.append(dict(method=meth,limit=lim,analyte=a,n=n,scheme=sc,k=len(te),cov=cov,
                           bias=100*np.median(te.e/T-1),hw_up=100*np.median(U/te.e-1)))
R=pd.DataFrame(rows); R.to_csv('work/ana/external_results.csv',index=False)
pd.set_option('display.width',220)
print(R.pivot_table(index=['method','analyte','limit','n'],columns='scheme',values='cov').round(2).to_string())
print(R[R.scheme=='native'].pivot_table(index=['method','analyte','limit'],columns='n',values='bias').round(1).to_string())
def summ(g): return pd.Series({'cells':len(g),'p90':100*(g['cov']>=0.9).mean(),'med':g['cov'].median(),'min':g['cov'].min()})
print(R.groupby(['scheme','method']).apply(summ,include_groups=False).round(3).to_string())
print(R[R.method.isin(['reflimR','RINet_B50'])].groupby('scheme').apply(summ,include_groups=False).round(3).to_string())
