# Symmetric vs signed (asymmetric) split-conformal recalibration, with and without shape-group Mondrian calibration,
# on the five person-disjoint splits (reflimR and RINet B=50). Also relative interval score (Winkler).
import pandas as pd, numpy as np
src=open('code/03_analysis/sens.py').read().split("rows=[]")[0]; exec(src)
A=0.05
def qhi(v,a): v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
def qlo(v,a): v=np.sort(v); k=int(np.floor((len(v)+1)*a)); return v[max(k,1)-1]
SHAPE={'TSH':'skew','Tg':'skew','GGT':'skew','TT4':'sym','TT3':'sym','HbA1c':'sym','FT4':'sym','Hb':'sym'}
rows=[]
for meth in ('reflimR','RINet_B50'):
  pL,pU,lLo,lHi,uLo,uHi=V[meth]
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['split','source','analyte','half','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['split','source','analyte','half','n','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['lest']=np.log(x.e); x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4)
    x['t']=(x.lest-np.log(x['T']))/x.s; x['sc']=x.t.abs(); x['shape']=x.analyte.map(SHAPE)
    S=x[x.source=='sim']
    for sp in range(1,6):
      Rr=x[(x.source=='nhanes')&(x.split==sp)]
      for th,ch in (('B','A'),('A','B')):
        for a in Rr.analyte.unique():
          for n in (300,500,1000):
            te=Rr[(Rr.half==th)&(Rr.analyte==a)&(Rr.n==n)]; cal=Rr[(Rr.half==ch)&(Rr.n==n)]
            base=pd.concat([S[S.n==n],cal[cal.analyte!=a]])            # scheme (d): analyte- and person-disjoint
            grp=base[base['shape']==SHAPE[a]]                             # shape-group Mondrian within (d)
            for cname,c in (('d',base),('d+shape',grp)):
              for kind in ('symmetric','signed'):
                if kind=='symmetric':
                  q=qhi(c.sc.values,A); L=np.exp(te.lest-q*te.s); U=np.exp(te.lest+q*te.s)
                else:
                  ql,qh=qlo(c.t.values,A/2),qhi(c.t.values,A/2); L=np.exp(te.lest-qh*te.s); U=np.exp(te.lest-ql*te.s)
                T=te['T']; cov=((L<=T+te.tol)&(T-te.tol<=U)).mean()
                IS=((U-L)+(2/A)*np.maximum(L-T,0)+(2/A)*np.maximum(T-U,0))/T
                rows.append(dict(split=sp,method=meth,limit=lim,half=th,analyte=a,n=n,calib=cname,kind=kind,cov=cov,
                                 hw_up=100*np.median(U/te.e-1),hw_down=100*np.median(1-L/te.e),relIS=np.median(IS),m=len(c)))
R=pd.DataFrame(rows); R.to_csv('work/ana/asym_results.csv',index=False)
def summ(g): return pd.Series({'p90':100*(g['cov']>=0.90).mean(),'p85':100*(g['cov']>=0.85).mean(),'med':g['cov'].median(),'min':g['cov'].min(),'relIS':g.relIS.median()})
P=R.groupby(['calib','kind','split']).apply(summ,include_groups=False)
print(P.groupby(['calib','kind']).agg(['mean','min','max']).round(3).to_string())
print(R.groupby(['calib','kind','method']).apply(summ,include_groups=False).round(3).to_string())
