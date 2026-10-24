// Port of reflimR 1.1.0 (Hoffmann, Klawitter, Trulson, Klawonn; GPL-3) core estimator with apply.rounding = FALSE.
function qnorm(p){ // Wichura AS241 (PPND16)
  if(p<=0) return -Infinity; if(p>=1) return Infinity;
  const q=p-0.5; let r,val;
  if(Math.abs(q)<=0.425){ r=0.180625-q*q;
    val=q*(((((((r*2509.0809287301226727+33430.575583588128105)*r+67265.770927008700853)*r+45921.953931549871457)*r+13731.693765509461125)*r+1971.5909503065514427)*r+133.14166789178437745)*r+3.387132872796366608)/
      (((((((r*5226.495278852545925+28729.085735721942674)*r+39307.89580009271061)*r+21213.794301586595867)*r+5394.1960214247511077)*r+687.1870074920579083)*r+42.313330701600911252)*r+1.0); return val; }
  r = q<0 ? p : 1-p; r=Math.sqrt(-Math.log(r));
  if(r<=5){ r-=1.6;
    val=(((((((r*7.7454501427834140764e-4+0.0227238449892691845833)*r+0.24178072517745061177)*r+1.27045825245236838258)*r+3.64784832476320460504)*r+5.7694972214606914055)*r+4.6303378461565452959)*r+1.42343711074968357734)/
      (((((((r*1.05075007164441684324e-9+5.475938084995344946e-4)*r+0.0151986665636164571966)*r+0.14810397642748007459)*r+0.68976733498510000455)*r+1.6763848301838038494)*r+2.05319162663775882187)*r+1.0);
  } else { r-=5;
    val=(((((((r*2.01033439929228813265e-7+2.71155556874348757815e-5)*r+0.0012426609473880784386)*r+0.026532189526576123093)*r+0.29656057182850489123)*r+1.7848265399172913358)*r+5.4637849111641143699)*r+6.6579046435011037772)/
      (((((((r*2.04426310338993978564e-15+1.4215117583164458887e-7)*r+1.8463183175100546818e-5)*r+7.868691311456132591e-4)*r+0.0148753612908506148525)*r+0.13692988092273580531)*r+0.59983220655588793769)*r+1.0); }
  return q<0 ? -val : val;
}
function quantile7(sorted,p){ const n=sorted.length, h=(n-1)*p, lo=Math.floor(h), hi=Math.min(lo+1,n-1); return sorted[lo]+(h-lo)*(sorted[hi]-sorted[lo]); }
const srt=a=>Float64Array.from(a).sort();
function bowley(x){ const s=srt(x), q1=quantile7(s,0.25), q2=quantile7(s,0.5), q3=quantile7(s,0.75); return (q1-2*q2+q3)/(q3-q1); }
function isLognormal(x){ const b1=bowley(x), b2=bowley(x.map(Math.log)); return b1<0 ? false : (b1-b2)>=0.05; }
function iboxplot(x,lognormal,percTrunc=2.5){
  let xx = lognormal ? x.map(Math.log) : x.slice(); const qt=percTrunc/100;
  const trunc=(v,qf)=>{ const s=srt(v), Q1=quantile7(s,0.25), Q2=quantile7(s,0.5), Q3=quantile7(s,0.75); const vr=Math.min(Q2-Q1,Q3-Q2);
    const lo=Q2-qf*vr, hi=Q2+qf*vr; return v.filter(e=>e>=lo&&e<=hi); };
  for(let i=1;i<=5;i++){ const tq=1-qt*i/5, al=(2*qt)*(i-1)/5; xx=trunc(xx, qnorm(tq)/qnorm(0.75-0.25*al)); }
  const qf=qnorm(1-qt)/qnorm(0.75-0.5*qt); let n0=1,n1=0;
  while(n0>n1){ n0=xx.length; xx=trunc(xx,qf); n1=xx.length; }
  return lognormal ? xx.map(Math.exp) : xx;
}
function truncatedQQ(xt,lognormal,percTrunc=2.5){
  const n=xt.length; let nq=100; if(n<nq) nq=n;
  const v= lognormal ? xt.map(Math.log) : xt; const s=srt(v);
  const xa=[],ya=[];
  for(let k=0;k<nq;k++){ const p1=percTrunc/100+(1-2*percTrunc/100)*k/(nq-1), p2=k/(nq-1); xa.push(qnorm(p1)); ya.push(quantile7(s,p2)); }
  const a0=Math.floor(0.05*nq), a1=Math.ceil(0.95*nq); // R 1-based inclusive indices
  let sx=0,sy=0,m=0; for(let i=a0;i<=a1;i++){ sx+=xa[i-1]; sy+=ya[i-1]; m++; }
  const mx=sx/m,my=sy/m; let sxy=0,sxx=0; for(let i=a0;i<=a1;i++){ const dx=xa[i-1]-mx; sxy+=dx*(ya[i-1]-my); sxx+=dx*dx; }
  const a=sxy/sxx, b=my-a*mx; let ll=b-1.96*a, ul=b+1.96*a;
  if(lognormal){ ll=Math.exp(ll); ul=Math.exp(ul); }
  if(ll<0) ll=0;
  return {mean:b, sd:a, lower:ll, upper:ul};
}
function confInt95(n,lower,upper,lognormal){
  let L=lower,U=upper; if(lognormal){ L=Math.log(L); U=Math.log(U); }
  const sg=(U-L)/3.92, dO=sg*5.81/(Math.sqrt(n)+0.66), dI=sg*7.26/(Math.sqrt(n)-5.58);
  let r=[L-dO, L+dI, U-dI, U+dO]; if(lognormal) r=r.map(Math.exp); return r;
}
function reflim(x){
  let xx=x.filter(v=>Number.isFinite(v)); const removed=xx.filter(v=>v<=0).length; xx=xx.filter(v=>v>0);
  const n=xx.length; if(n<40) return {ok:false, reason:'n<40', n};
  const lognormal=isLognormal(xx); const tr=iboxplot(xx,lognormal);
  if(tr.length<40) return {ok:false, reason:'n<40 after truncation', n, nTrunc:tr.length};
  const qq=truncatedQQ(tr,lognormal);
  if(!(qq.upper>qq.lower)) return {ok:false, reason:'invalid limits', n};
  const ci=confInt95(n,qq.lower,qq.upper,lognormal);
  return {ok:true, n, nTrunc:tr.length, lognormal, removed, lower:qq.lower, upper:qq.upper, ciLower:[ci[0],ci[1]], ciUpper:[ci[2],ci[3]]};
}
if(typeof module!=='undefined') module.exports={reflim,qnorm};
