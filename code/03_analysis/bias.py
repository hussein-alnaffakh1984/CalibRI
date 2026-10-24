import pandas as pd, numpy as np
exec(open('code/03_analysis/disjoint.py').read().split("rows=[]")[0])
B=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  if meth.startswith('refineR'): meth='refineR'
  for lim,p,t in (('LRL',pL,'trueLL'),('URL',pU,'trueUL')):
    x=d[(d.source=='nhanes')&d[p].notna()&(d[p]>0)]
    g=x.assign(rb=(x[p]/x[t]-1)*100).groupby(['analyte','n_nominal','half']).rb.median().groupby(['analyte','n_nominal']).mean()
    for (a,n),v in g.items(): B.append(dict(method=meth,limit=lim,analyte=a,n=n,bias=v))
B=pd.DataFrame(B); B.to_csv('work/ana/bias.csv',index=False)
P=B.pivot_table(index=['analyte','limit'],columns=['method','n'],values='bias').round(1)
print(P.to_string())
