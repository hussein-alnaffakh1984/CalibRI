import pandas as pd, numpy as np
src=open('code/03_analysis/sens.py').read().split("rows=[]")[0]; exec(src)
rows=[]
for meth in ('reflimR','RINet_B50'):
  pL,pU,lLo,lHi,uLo,uHi=V[meth]
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['split','source','analyte','half','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['split','source','analyte','half','n','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(x.lest-np.log(x['T']))/x.s
    S=x[x.source=='sim']
    for sp in range(1,6):
      Rr=x[(x.source=='nhanes')&(x.split==sp)]
      for th,ch in (('B','A'),('A','B')):
        for a in Rr.analyte.unique():
          for n in (300,500,1000):
            te=Rr[(Rr.half==th)&(Rr.analyte==a)&(Rr.n==n)]; cal=Rr[(Rr.half==ch)&(Rr.n==n)]
            for sc,c in (('native',None),('a',S[S.n==n]),('c',pd.concat([S[S.n==n],cal])),('d',pd.concat([S[S.n==n],cal[cal.analyte!=a]])),('d real-only',cal[cal.analyte!=a]),('c real-only',cal)):
              if c is None: L,U=te.lo,te.hi
              else: q=qhat(c.sc.values); L,U=np.exp(te.lest-q*te.s),np.exp(te.lest+q*te.s)
              T=te['T']; IS=((U-L)+(2/0.05)*np.maximum(L-T,0)+(2/0.05)*np.maximum(T-U,0))/T
              rows.append(dict(split=sp,method=meth,scheme=sc,cov=((L<=T+te.tol)&(T-te.tol<=U)).mean(),IS=np.median(IS),ISmean=np.mean(IS)))
R=pd.DataFrame(rows); R.to_csv('work/ana/is_all.csv',index=False)
P=R.groupby(['scheme','split']).apply(lambda g: pd.Series({'p90':100*(g['cov']>=0.9).mean(),'med':g['cov'].median(),'min':g['cov'].min(),'IS':g.IS.median(),'ISmean':g.ISmean.mean()}),include_groups=False)
print(P.groupby('scheme').agg(['mean','min','max']).round(3).to_string())
