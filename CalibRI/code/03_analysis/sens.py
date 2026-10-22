# Sensitivity analyses: (1) five person-disjoint splits; (2) RINet bootstrap with 50 vs 1,000 resamples.
import pandas as pd, numpy as np, sys
TOL={'TT4':0.05,'TT3':0.5,'HbA1c':0.05}
def qhat(v,a=0.05):
    v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
# ---- simulated calibration pool (n = 300/500/1000)
m1=pd.read_csv('work/ribundle/meta.csv'); m1=m1[(m1.source=='sim')&m1.n_nominal.isin([300,500,1000])]
f1=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv']])
rs=pd.read_csv('work/res_b4/rinet_sim.csv')
sim=m1.merge(f1[['ds_id','rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi']],on='ds_id').merge(rs,on='ds_id')
sim['half']='S'; sim['split']=0
# ---- test data: split 1 (bundle 3) and splits 2-5 (bundle 4)
m3=pd.read_csv('work/ribundle3/meta.csv'); m3['split']=1
f3=pd.read_csv('work/res_b3/res_fast_shard01.csv')[['ds_id','rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi']]
r3=pd.read_csv('work/res_b4/rinet_b3.csv')
s1=m3.merge(f3,on='ds_id').merge(r3,on='ds_id')
m4=pd.read_csv('work/ribundle4/meta.csv'); f4=pd.read_csv('work/res_b4/reflim_b4.csv'); r4=pd.read_csv('work/res_b4/rinet_b4.csv')
s4=m4.merge(f4,on='ds_id').merge(r4,on='ds_id')
real=pd.concat([s1,s4],ignore_index=True)
print('datasets: sim',len(sim),'real',len(real),'per split',real.groupby('split').size().to_dict())
d=pd.concat([sim,real],ignore_index=True)
d['tol']=np.where(d.source=='nhanes',d.analyte.map(TOL).fillna(0),0.0)
V={'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet_B50':('rn_LL','rn_UL','B50_LL_lo','B50_LL_hi','B50_UL_lo','B50_UL_hi'),
   'RINet_B1000':('rn_LL','rn_UL','B1000_LL_lo','B1000_LL_hi','B1000_UL_lo','B1000_UL_hi')}
rows=[]; widths=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in V.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['split','source','analyte','half','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['split','source','analyte','half','n','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(x.lest-np.log(x['T']))/x.s
    S=x[x.source=='sim']
    for sp in sorted(x[x.source=='nhanes'].split.unique()):
      R_=x[(x.source=='nhanes')&(x.split==sp)]
      for th,ch in (('B','A'),('A','B')):
        for a in R_.analyte.unique():
          for n in (300,500,1000):
            te=R_[(R_.half==th)&(R_.analyte==a)&(R_.n==n)]; cal=R_[(R_.half==ch)&(R_.n==n)]
            for scheme in ('native','sim-only','sim+otherHalf(all)','sim+otherHalf(other)'):
              if scheme=='native': L,U=te.lo,te.hi; hwv=np.sqrt(te.hi/te.lo)-1
              else:
                c=S[S.n==n]
                if scheme=='sim+otherHalf(all)': c=pd.concat([c,cal])
                if scheme=='sim+otherHalf(other)': c=pd.concat([c,cal[cal.analyte!=a]])
                q=qhat(c.sc.values); L,U=np.exp(te.lest-q*te.s),np.exp(te.lest+q*te.s); hwv=np.exp(q*te.s)-1
              cov=((L<=te['T']+te.tol)&(te['T']-te.tol<=U)).mean()
              rows.append(dict(split=sp,method=meth,limit=lim,half=th,analyte=a,n=n,scheme=scheme,cov=cov,hw=100*np.median(hwv)))
R=pd.DataFrame(rows); R.to_csv('work/ana/sens_results.csv',index=False)
def summ(g):
    return pd.Series({'cells':len(g),'pct_ge90':100*(g['cov']>=0.90).mean(),'pct_ge85':100*(g['cov']>=0.85).mean(),'median':g['cov'].median(),'min':g['cov'].min()})
# (1) Table-2 analogue per split: reflimR + RINet (B=50, seeded per dataset)
A=R[R.method.isin(['reflimR','RINet_B50'])].groupby(['scheme','split']).apply(summ,include_groups=False)
print('\n== reflimR + RINet(B50): per split'); print(A.round(3).to_string())
print('\n== across 5 splits: mean [min, max]')
for sc,g in A.groupby(level=0):
    print(f"{sc:24s}", '  '.join(f"{c}={g[c].mean():.2f} [{g[c].min():.2f},{g[c].max():.2f}]" for c in ['pct_ge90','pct_ge85','median','min']))
# (2) RINet B50 vs B1000, pooled over splits
Bc=R[R.method.str.startswith('RINet')].groupby(['method','scheme']).apply(summ,include_groups=False)
print('\n== RINet B=50 vs B=1000 (all 5 splits pooled)'); print(Bc.round(3).to_string())
Wd=R[R.method.str.startswith('RINet')&(R.limit=='URL')].groupby(['method','scheme','analyte','n']).hw.mean().unstack(['method']).round(1)
print('\n== RINet URL median upward half-width %, native and sim+otherHalf(all)'); print(Wd.loc[['native','sim+otherHalf(all)']].to_string())
