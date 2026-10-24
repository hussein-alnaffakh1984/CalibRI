import pandas as pd, numpy as np
exec(open('code/03_analysis/disjoint.py').read().split("rows=[]")[0])
W=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('URL',pU,uLo,uHi,'trueUL'),('LRL',pL,lLo,lHi,'trueLL')):
    x=d[['source','analyte','half','n_nominal','tol',p,lo,hi,t]].copy(); x.columns=['source','analyte','half','n','tol','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(np.log(x.e)-np.log(x['T']))/x.s
    for th,ch in (('B','A'),('A','B')):
      for a in ['TSH','Tg','TT4','TT3','HbA1c']:
        for n in (300,500,1000):
          c=x[((x.source=='sim')|((x.source=='nhanes')&(x.half==ch)&(x.analyte!=a)))&(x.n==n)]; te=x[(x.source=='nhanes')&(x.half==th)&(x.analyte==a)&(x.n==n)]
          q=qhat(c.sc.values); W.append(dict(method=meth,limit=lim,analyte=a,n=n,half=th,hw=np.median(np.exp(q*te.s)-1)*100,native_hw=np.median(np.sqrt(te.hi/te.lo)-1)*100))
W=pd.DataFrame(W).groupby(['method','limit','analyte','n'])[['hw','native_hw']].mean().reset_index()
W['method']=W.method.replace({'refineR(UM-scale)':'refineR'}); W.to_csv('work/ana/widths_d.csv',index=False)
print(W[W.limit=='URL'].pivot_table(index=['analyte','n'],columns='method',values='hw').round(0).to_string())
print(W[W.limit=='LRL'].pivot_table(index=['analyte','n'],columns='method',values='hw').round(0).to_string())
