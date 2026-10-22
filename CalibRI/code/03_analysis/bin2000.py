# Check of the n = 2,000 calibration bin on person-disjoint data (split 1), reflimR and RINet (B = 50).
import pandas as pd, numpy as np
TOL={'TT4':0.05,'TT3':0.5,'HbA1c':0.05}
def qhat(v,a=0.05):
    v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
m1=pd.read_csv('work/ribundle/meta.csv'); m1=m1[(m1.source=='sim')&(m1.n_nominal==2000)]
f1=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv']])
r1=pd.read_csv('work/res_C/res_rinet.csv')   # B=50 single-stream run (sim n=2000 not in seeded re-run)
sim=m1.merge(f1,on='ds_id').merge(r1,on='ds_id'); sim['src']='sim'
sim['B50_LL_lo']=sim.rn_CI95_LL_lo; sim['B50_LL_hi']=sim.rn_CI95_LL_hi; sim['B50_UL_lo']=sim.rn_CI95_UL_lo; sim['B50_UL_hi']=sim.rn_CI95_UL_hi
m7=pd.read_csv('work/ribundle7/meta.csv'); nh=m7.merge(pd.read_csv('work/res_b7/reflim_b7.csv'),on='ds_id').merge(pd.read_csv('work/res_b7/rinet_b7.csv'),on='ds_id'); nh['src']='nhanes'
d=pd.concat([sim,nh],ignore_index=True); d['tol']=np.where(d.src=='nhanes',d.analyte.map(TOL).fillna(0),0.0)
V={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','B50_LL_lo','B50_LL_hi','B50_UL_lo','B50_UL_hi')}
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in V.items():
  for lim,p,lo,hi,t in (('LRL',pL,lLo,lHi,'trueLL'),('URL',pU,uLo,uHi,'trueUL')):
    x=d[['src','analyte','half','tol',p,lo,hi,t]].copy(); x.columns=['src','analyte','half','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(x.lest-np.log(x['T']))/x.s
    S=x[x.src=='sim']
    for th,ch in (('B','A'),('A','B')):
      for a in x[x.src=='nhanes'].analyte.unique():
        te=x[(x.src=='nhanes')&(x.half==th)&(x.analyte==a)]; cal=x[(x.src=='nhanes')&(x.half==ch)]
        for sc,c in (('native',None),('c',pd.concat([S,cal])),('d',pd.concat([S,cal[cal.analyte!=a]]))):
          if c is None: L,U=te.lo,te.hi
          else: q=qhat(c.sc.values); L,U=np.exp(te.lest-q*te.s),np.exp(te.lest+q*te.s)
          T=te['T']; rows.append(dict(method=meth,limit=lim,analyte=a,half=th,scheme=sc,cov=((L<=T+te.tol)&(T-te.tol<=U)).mean()))
R=pd.DataFrame(rows); R.to_csv('work/ana/bin2000_results.csv',index=False)
print(R.groupby('scheme').apply(lambda g: pd.Series({'cells':len(g),'p90':100*(g['cov']>=0.9).mean(),'med':g['cov'].median(),'min':g['cov'].min()}),include_groups=False).round(3))
print(R[R.analyte=='TSH'].pivot_table(index=['method','limit','half'],columns='scheme',values='cov').round(2))
