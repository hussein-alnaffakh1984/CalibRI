// RINet 1D (LeBien, Velev & Roche-Lima; rinet R package 0.1.1, MIT licence) — forward pass of the bundled CNN in plain JavaScript.
const RINET=(function(){
  let W=null, CENTER=null, SCALE=null;
  function load(spec){
    const bin = typeof atob==='function' ? Uint8Array.from(atob(spec.b64),c=>c.charCodeAt(0)) : Uint8Array.from(Buffer.from(spec.b64,'base64'));
    const all=new Float32Array(bin.buffer, bin.byteOffset, bin.byteLength/4); let o=0; W=[];
    for(const sh of spec.shapes){ const n=sh.reduce((a,b)=>a*b,1); W.push(all.subarray(o,o+n)); o+=n; }
    CENTER=spec.center; SCALE=spec.scale; }
  function conv(x,L,cin,k,b,cout){ // same padding, stride 1, kernel 5, relu; x[t*cin+c], k[(j*cin+c)*cout+o]
    const y=new Float32Array(L*cout);
    for(let t=0;t<L;t++) for(let o=0;o<cout;o++){ let s=b[o];
      for(let j=0;j<5;j++){ const tt=t+j-2; if(tt<0||tt>=L) continue; const xb=tt*cin, kb=j*cin;
        for(let c=0;c<cin;c++) s+=x[xb+c]*k[(kb+c)*cout+o]; }
      y[t*cout+o]= s>0?s:0; }
    return y; }
  function pool(x,L,C){ const Lo=Math.ceil(L/2), y=new Float32Array(Lo*C);
    for(let t=0;t<Lo;t++) for(let c=0;c<C;c++){ const a=x[(2*t)*C+c]; const b= 2*t+1<L ? x[(2*t+1)*C+c] : -Infinity; y[t*C+c]= a>b?a:b; }
    return [y,Lo]; }
  function dense(x,w,b,nin,nout,relu){ const y=new Float32Array(nout);
    for(let o=0;o<nout;o++){ let s=b[o]; for(let i=0;i<nin;i++) s+=x[i]*w[i*nout+o]; y[o]= relu ? (s>0?s:0) : s; } return y; }
  function predict(h){ // h: 100 features
    let x=Float32Array.from(h), L=100;
    x=conv(x,L,1,W[0],W[1],32); [x,L]=pool(x,L,32);
    x=conv(x,L,32,W[2],W[3],64); [x,L]=pool(x,L,64);
    x=conv(x,L,64,W[4],W[5],64); [x,L]=pool(x,L,64);
    x=dense(x,W[6],W[7],832,32,true); x=dense(x,W[8],W[9],32,3,false);
    return [x[0]*SCALE[0]+CENTER[0], x[1]*SCALE[1]+CENTER[1], x[2]*SCALE[2]+CENTER[2]]; }
  function features(lx){ const n=lx.length; let m=0; for(const v of lx) m+=v; m/=n; let ss=0; for(const v of lx) ss+=(v-m)*(v-m); const s=Math.sqrt(ss/n);
    const h=new Array(100).fill(0);
    for(const v of lx){ const z=(v-m)/s; if(z<-4||z>4) continue; let i=Math.floor((z+4)/0.08); if(i>99) i=99;
      // match numpy edges exactly: edges are linspace(-4,4,101)
      const lo=-4+8*i/100, hi=-4+8*(i+1)/100; if(z<lo&&i>0) i--; else if(z>=hi&&i<99) i++;
      h[i]++; }
    let mn=Infinity,mx=-Infinity; for(const c of h){ if(c<mn) mn=c; if(c>mx) mx=c; }
    return {h:h.map(c=>(c-mn)/(mx-mn)), m, s}; }
  const Z=1.959963984540054;
  function point(lx){ const f=features(lx); const p=predict(f.h); const mu=p[0]*f.s+f.m, sd=p[1]*f.s; return [Math.exp(mu-Z*sd), Math.exp(mu+Z*sd)]; }
  function q7(arr,p){ const s=Float64Array.from(arr).sort(); const h=(s.length-1)*p, lo=Math.floor(h), hi=Math.min(lo+1,s.length-1); return s[lo]+(h-lo)*(s[hi]-s[lo]); }
  function rng32(a){ return function(){ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
  function estimate(x,B=50,seed=20261016){
    let v=x.filter(Number.isFinite); const pos=v.filter(e=>e>0); if(pos.length<40) return {ok:false, reason:'n<40', n:v.length};
    const nonpos=v.length-pos.length; if(nonpos){ const mn=Math.min(...pos)/2; v=v.map(e=>e>0?e:mn); }
    const lx=v.map(Math.log); const [LL,UL]=point(lx);
    const r=rng32(seed), n=lx.length, bl=[], bu=[], s=new Float64Array(n);
    for(let b=0;b<B;b++){ for(let i=0;i<n;i++) s[i]=lx[Math.floor(r()*n)]; const [l,u]=point(s); bl.push(l); bu.push(u); }
    return {ok:true, n:v.length, nonpos, lower:LL, upper:UL, ciLower:[q7(bl,0.025),q7(bl,0.975)], ciUpper:[q7(bu,0.025),q7(bu,0.975)]}; }
  return {load, predict, features, point, estimate};
})();
if(typeof module!=='undefined') module.exports=RINET;
