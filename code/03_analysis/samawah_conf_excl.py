import pandas as pd, numpy as np
# calibration pool: sim (bundle1) + all real NHANES arms except coarse FT4 (bundles 1,2,3)
m=pd.concat([pd.read_csv('work/ribundle/meta.csv'),pd.read_csv('work/ribundle2/meta.csv'),pd.read_csv('work/ribundle3/meta.csv')])
f=pd.concat([pd.read_csv(p) for p in ['work/res_A1/res_fast_shard01.csv','work/res_A2/res_fast_shard02.csv','work/res_b2/res_fast_shard01.csv','work/res_b3/res_fast_shard01.csv']])
r=pd.concat([pd.read_csv(p) for p in ['work/res_C/res_rinet.csv','work/res_b2/res_rinet.csv','work/res_b3/res_rinet.csv']])
cal=m.merge(f,on='ds_id').merge(r,on='ds_id'); cal=cal[~((cal.source=='nhanes')&(cal.analyte=='FT4'))]
sm=pd.read_csv('work/samawah_excl/meta.csv').merge(pd.read_csv('work/samawah_excl/out/res_fast_shard01.csv'),on='ds_id').merge(pd.read_csv('work/samawah_excl/out/res_rinet.csv'),on='ds_id')
M={'refineR':('rf_LL','rf_UL','rf_UM95_LL_lo','rf_UM95_LL_hi','rf_UM95_UL_lo','rf_UM95_UL_hi'),
   'reflimR':('rl_LL','rl_UL','rl_CI95_LL_lo','rl_CI95_LL_hi','rl_CI95_UL_lo','rl_CI95_UL_hi'),
   'RINet':('rn_LL','rn_UL','rn_CI95_LL_lo','rn_CI95_LL_hi','rn_CI95_UL_lo','rn_CI95_UL_hi')}
BINS=np.array([300,500,1000,2000])
def qhat(v,a=0.05):
    v=np.sort(v); k=int(np.ceil((len(v)+1)*(1-a))); return v[min(k,len(v))-1]
out=[]
for meth,(pL,pU,lLo,lHi,uLo,uHi) in M.items():
  for lim,p,lo,hi,t in (('LRL',pL,lLo,lHi,'trueLL'),('URL',pU,uLo,uHi,'trueUL')):
    c=cal[['n_nominal',p,lo,hi,t]].dropna(); c=c[(c[p]>0)&(c[lo]>0)&(c[hi]>0)]
    s=np.maximum((np.log(c[hi])-np.log(c[lo]))/(2*1.96),1e-4); c=c.assign(sc=np.abs(np.log(c[p])-np.log(c[t]))/s)
    for _,row in sm.iterrows():
        b=BINS[np.argmin(np.abs(np.log(BINS)-np.log(row.n)))]
        q=qhat(c[c.n_nominal==b].sc.values); ss=max((np.log(row[hi])-np.log(row[lo]))/(2*1.96),1e-4) if row[lo]>0 else np.nan
        out.append(dict(method=meth,limit=lim,sex=row.sex,period=row.period,n=row.n,bin=b,est=row[p],native_lo=row[lo],native_hi=row[hi],q=q,
                        cal_lo=row[p]*np.exp(-q*ss),cal_hi=row[p]*np.exp(q*ss)))
O=pd.DataFrame(out); O.to_csv('work/samawah_excl/samawah_excl_results.csv',index=False)
pd.set_option('display.width',250)
print(O.round(3).to_string(index=False))
