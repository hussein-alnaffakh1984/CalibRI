import pandas as pd, numpy as np
# simulated calibration pool (bundle 1, sim arms)
m1=pd.read_csv('work/ribundle/meta.csv'); m1=m1[m1.source=='sim']; m1['half']='S'
f1=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv']]); r1=pd.read_csv('work/res_C/res_rinet.csv')
s=m1.merge(f1,on='ds_id').merge(r1,on='ds_id')
m3=pd.read_csv('work/ribundle3/meta.csv'); f3=pd.read_csv('work/res_b3/res_fast_shard01.csv'); r3=pd.read_csv('work/res_b3/res_rinet.csv')
h=m3.merge(f3,on='ds_id').merge(r3,on='ds_id')
d=pd.concat([s,h],ignore_index=True); d=d[d.n_nominal.isin([300,500,1000])]
TOL={'TT4':0.05,'TT3':0.5,'HbA1c':0.05}
d['tol']=np.where(d.source=='nhanes',d.analyte.map(TOL).fillna(0),0.0)
M={'refineR(UM-scale)':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi'),
   'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','rn_CI95_LL_lo','rn_CI95_LL_hi','rn_CI95_UL_lo','rn_CI95_UL_hi')}
def qhat(v,a=0.05):
    v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['source','analyte','half','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['source','analyte','half','n','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(x.lest-np.log(x['T']))/x.s
    real=x[x.source=='nhanes']
    for test_half,cal_half in (('B','A'),('A','B')):
      for a in real.analyte.unique():
        te_all=real[(real.half==test_half)&(real.analyte==a)]
        for scheme in ('native','sim-only','sim+otherHalf(all analytes)','sim+otherHalf(other analytes)'):
          for n in (300,500,1000):
            te=te_all[te_all.n==n]
            for tm in ('tol','notol'):
              tt=te.tol if tm=='tol' else 0
              if scheme=='native':
                cov=((te.lo<=te['T']+tt)&(te['T']-tt<=te.hi)).mean()
              else:
                if scheme=='sim-only': c=x[(x.source=='sim')&(x.n==n)]
                elif scheme=='sim+otherHalf(all analytes)': c=x[((x.source=='sim')|((x.source=='nhanes')&(x.half==cal_half)))&(x.n==n)]
                else: c=x[((x.source=='sim')|((x.source=='nhanes')&(x.half==cal_half)&(x.analyte!=a)))&(x.n==n)]
                q=qhat(c.sc.values); hw=q*te.s
                cov=((np.exp(te.lest-hw)<=te['T']+tt)&(te['T']-tt<=np.exp(te.lest+hw))).mean()
              rows.append(dict(method=meth,limit=lim,test_half=test_half,analyte=a,scheme=scheme,n=n,tol=tm,cov=cov))
R=pd.DataFrame(rows); R.to_csv('work/ana/disjoint_results.csv',index=False)
pd.set_option('display.width',250)
def summ(sub,label):
    g=sub.groupby(['scheme','tol'])['cov'].agg(cells='size',ge90=lambda v:(v>=.9).mean(),ge85=lambda v:(v>=.85).mean(),mn='min',med='median').round(3)
    print('\n#####',label); print(g.to_string())
summ(R[R.method.isin(['reflimR','RINet'])],'GENUINE CIs (reflimR + RINet): 5 analytes x 2 halves x 3 n x 2 limits x 2 methods = 120 cells')
summ(R[R.method=='refineR(UM-scale)'],'refineR with UM as scale (UM not an estimation CI)')
w=R[(R.scheme=='sim+otherHalf(other analytes)')&(R.tol=='tol')&R.method.isin(['reflimR','RINet'])]
print('\nworst cells, strictest scheme (genuine CIs, tol):'); print(w.sort_values('cov').head(10).to_string(index=False))
