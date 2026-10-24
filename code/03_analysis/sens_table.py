import pandas as pd, numpy as np, json
from decimal import Decimal, ROUND_HALF_UP
def r0(x): return str(Decimal(str(x)).quantize(Decimal('1'),ROUND_HALF_UP))
def r2(x): return str(Decimal(str(round(x,6))).quantize(Decimal('0.01'),ROUND_HALF_UP))
R=pd.read_csv('sens_results.csv')
def summ(g): return pd.Series({'p90':100*(g['cov']>=0.90).mean(),'p85':100*(g['cov']>=0.85).mean(),'med':g['cov'].median(),'min':g['cov'].min()})
lab={'native':'Native','sim-only':'Conformal: simulation only','sim+otherHalf(all)':'Conformal: simulation + other-half real, all analytes','sim+otherHalf(other)':'Conformal: simulation + other-half real, other analytes'}
order=['native','sim-only','sim+otherHalf(all)','sim+otherHalf(other)']
A=R[R.method.isin(['reflimR','RINet_B50'])].groupby(['scheme','split']).apply(summ,include_groups=False)
rowsA=[]
for sc in order:
    g=A.loc[sc]
    f=lambda c,fn: f"{fn(g[c].mean())} ({fn(g[c].min())}–{fn(g[c].max())})"
    rowsA.append([lab[sc],f('p90',r0)+'%' if False else f"{r0(g.p90.mean())}% ({r0(g.p90.min())}–{r0(g.p90.max())}%)",f"{r0(g.p85.mean())}% ({r0(g.p85.min())}–{r0(g.p85.max())}%)",f"{r2(g.med.mean())} ({r2(g.med.min())}–{r2(g.med.max())})",f"{r2(g['min'].mean())} ({r2(g['min'].min())}–{r2(g['min'].max())})"])
B=R[R.method.str.startswith('RINet')].groupby(['scheme','method']).apply(summ,include_groups=False)
rowsB=[]
for sc in order:
    a,b=B.loc[(sc,'RINet_B50')],B.loc[(sc,'RINet_B1000')]
    rowsB.append([lab[sc],f"{r0(a.p90)}% / {r0(b.p90)}%",f"{r0(a.p85)}% / {r0(b.p85)}%",f"{r2(a.med)} / {r2(b.med)}",f"{r2(a['min'])} / {r2(b['min'])}"])
json.dump({'A':rowsA,'B':rowsB},open('sens_table.json','w'),ensure_ascii=False,indent=1)
for r in rowsA+rowsB: print(' | '.join(r))
