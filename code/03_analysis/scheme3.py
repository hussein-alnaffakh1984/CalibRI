import pandas as pd, numpy as np
exec(open('code/03_analysis/full.py').read().split("nat=[]; conf=[]")[0])   # reuse data loading only
THY=['nhanes:TSH','nhanes:TT4','nhanes:TT3','nhanes:Tg']; REAL=THY+['nhanes:HbA1c']
def qhat(s,a=0.05):
    s=np.sort(s); m=len(s); k=int(np.ceil((m+1)*(1-a))); return s[min(k,m)-1]
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['arm','source','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['arm','source','n','tol','e','lo','hi','T']
    y=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy(); y=y[y.arm!='nhanes:FT4']
    y['lest']=np.log(y.e); y['s']=np.maximum((np.log(y.hi)-np.log(y.lo))/(2*1.96),1e-4); y['nrm']=np.abs(y.lest-np.log(y['T']))/y.s
    for a in REAL:
      other = ['nhanes:HbA1c'] if a in THY else THY      # different survey component / mostly different people
      for tolmode in ('tol','notol'):
        for n in sorted(y.n.unique()):
          c=y[((y.source=='sim')|y.arm.isin(other))&(y.n==n)]; te=y[(y.arm==a)&(y.n==n)]
          q=qhat(c.nrm.values); h=q*te.s; tt=te.tol if tolmode=='tol' else 0
          rows.append(dict(method=meth,limit=lim,arm=a,n=n,tolmode=tolmode,cov=((np.exp(te.lest-h)<=te['T']+tt)&(te['T']-tt<=np.exp(te.lest+h))).mean()))
R=pd.DataFrame(rows)
for tm in ('tol','notol'):
  r=R[R.tolmode==tm]; print(tm,'cells>=0.90',round((r['cov']>=.9).mean(),3),'>=0.85',round((r['cov']>=.85).mean(),3),'min',round(r['cov'].min(),2))
  print(r[r['cov']<0.85].sort_values('cov')[['method','limit','arm','n','cov']].head(12).to_string(index=False))
