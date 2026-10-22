import pandas as pd, numpy as np, json
from decimal import Decimal, ROUND_HALF_UP
def r(x,q): return str(Decimal(str(round(float(x),9))).quantize(Decimal(q),ROUND_HALF_UP)).replace('-','−')
src=open('code/03_analysis/disjoint.py').read().split("rows=[]")[0]; exec(src)
rows=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  name={'refineR(UM-scale)':'refineR (UM)','reflimR':'reflimR','RINet':'RINet'}[meth]
  for lim,p,lo,hi,t in (('LRL',pL,lLo,lHi,'trueLL'),('URL',pU,uLo,uHi,'trueUL')):
    x=d[['source','analyte','half','n_nominal',p,lo,hi,t]].copy(); x.columns=['source','analyte','half','n','e','lo','hi','T']
    x=x[x[['e','lo','hi']].notna().all(axis=1)&(x.e>0)&(x.lo>0)&(x.hi>0)].copy()
    x['s']=np.maximum((np.log(x.hi)-np.log(x.lo))/(2*1.96),1e-4); x['sc']=np.abs(np.log(x.e)-np.log(x['T']))/x.s
    for n in (300,500,1000):
      ms=[];mr=[];qs=[]
      for th,ch in (('B','A'),('A','B')):
        for a in ['TSH','Tg','TT4','TT3','HbA1c']:
          S=x[(x.source=='sim')&(x.n==n)]; Rr=x[(x.source=='nhanes')&(x.half==ch)&(x.n==n)&(x.analyte!=a)]
          ms.append(len(S)); mr.append(len(Rr)); qs.append(qhat(pd.concat([S,Rr]).sc.values))
      rows.append([name,lim,f'{n:,}',f'{min(ms):,}',f'{min(mr):,}–{max(mr):,}',f'{r(min(qs),"0.01")}–{r(max(qs),"0.01")}'])
E=pd.read_csv('work/ana/external_results.csv')
ext=[]
for m,lab in (('reflimR','reflimR'),('RINet_B50','RINet (B = 50)'),('RINet_B1000','RINet (B = 1,000)'),('refineR_UM','refineR (UM)')):
  for a in ('TSH','FT4'):
    for lim in ('LRL','URL'):
      g=E[(E.method==m)&(E.analyte==a)&(E.limit==lim)]
      nat=g[g.scheme=='native'].sort_values('n'); dd=g[g.scheme=='sim+NHANES(other)'].sort_values('n')
      ext.append([lab,a,lim,' / '.join(r(v,'1') for v in nat.bias),' / '.join(r(v,'0.01') for v in nat['cov']),' / '.join(r(v,'0.01') for v in dd['cov']),' / '.join(r(v,'1') for v in dd.hw_up)])
json.dump({'s8':ext,'s9':rows},open('work/ana/supp_extra.json','w'),ensure_ascii=False)
for x_ in rows[:4]: print(x_)
for x_ in ext[:4]: print(x_)
