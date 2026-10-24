const BINS=[300,500,1000,2000];
function nearestBin(n){ let b=BINS[0],d=Infinity; for(const x of BINS){ const e=Math.abs(Math.log(x)-Math.log(n)); if(e<d){d=e;b=x;} } return b; }
function calibrate(est,lo,hi,q){ if(!(lo>0&&hi>0&&est>0)) return null; const s=Math.max((Math.log(hi)-Math.log(lo))/(2*1.96),1e-4); return [est*Math.exp(-q*s), est*Math.exp(q*s)]; }
if(typeof module!=='undefined') module.exports={nearestBin,calibrate,BINS};
