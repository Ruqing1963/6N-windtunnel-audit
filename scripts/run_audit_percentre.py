import numpy as np, json, time, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); os.makedirs(DATA, exist_ok=True)
t0=time.time()
LIMIT=200_000_000; NMAX=LIMIT//6
sieve=np.ones(LIMIT+3,dtype=bool); sieve[:2]=False
for p in range(2,int((LIMIT+2)**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
is_twin=np.zeros(NMAX+1,dtype=bool)
for a in range(1,NMAX+1,10_000_000):
    b=min(a+10_000_000,NMAX+1); idx=np.arange(a,b,dtype=np.int64)
    is_twin[a:b]=sieve[6*idx-1]&sieve[6*idx+1]
ntwin=int(is_twin.sum()); p0=ntwin/NMAX

omega=np.zeros(NMAX+1,dtype=np.int8)
M  =np.ones(NMAX+1,dtype=np.float64)   # wind-tunnel object: prod q/(q-2)
M2 =np.ones(NMAX+1,dtype=np.float64)   # correct per-centre core: prod (q-1)/(q-3)
primes=np.nonzero(sieve[:NMAX+1])[0]; primes=primes[primes>=5]
for q in primes:
    omega[q::q]+=1; M[q::q]*= q/(q-2); M2[q::q]*= (q-1)/(q-3)

# global constant C = prod_{q>3} q(q-3)/((q-1)(q-2)), over all relevant primes (converges)
qq=primes.astype(np.float64)
C=np.prod( qq*(qq-3)/((qq-1)*(qq-2)) )
print("global suppression C = %.6f  (primes up to %d)"%(C, primes[-1]))
M2 *= C   # full per-centre model:  E[twin|N]/p0  ~=  C * prod_{q|N}(q-1)/(q-3)

# ---- TEST: does the CORRECT per-centre model M2 predict measured enrichment 1:1? ----
edges=np.array([0.6,0.8,0.95,1.1,1.3,1.6,2.0,2.6,4.0])
mc=M2[1:]; tw=is_twin[1:]; wm=M[1:]
print("\n  bin (model M2)      n        <M2>    <M_wind>   meas_enrich   meas/M2   meas/M_wind")
binrows=[]
for i in range(len(edges)-1):
    sel=(mc>=edges[i])&(mc<edges[i+1]); n=int(sel.sum())
    if n<2000: continue
    meas=float(tw[sel].mean()/p0); m2=float(mc[sel].mean()); mw=float(wm[sel].mean())
    binrows.append(dict(lo=float(edges[i]),hi=float(edges[i+1]),n=n,M2=m2,Mwind=mw,meas=meas,
                        r_M2=meas/m2,r_wind=meas/mw))
    print("  [%.2f,%.2f)  %9d   %.4f   %.4f    %.4f      %.4f   %.4f"%(
        edges[i],edges[i+1],n,m2,mw,meas,meas/m2,meas/mw))

# ---- primorial centre: wind-tunnel x7.03 vs correct per-centre value ----
P=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]
wind=1.0; core=1.0
for q in P: wind*=q/(q-2); core*=(q-1)/(q-3)
percentre=C*core
print("\nPRIMORIAL Omega_15:")
print("  wind-tunnel  prod q/(q-2)            = %.4f   (the x7.03 figure)"%wind)
print("  correct per-centre  C*prod(q-1)/(q-3)= %.4f"%percentre)
print("  ratio per-centre / wind-tunnel       = %.4f  (%.2f%% gap)"%(percentre/wind,100*(percentre/wind-1)))

# ---- per-omega measured vs correct model <M2> ----
om=[]
for k in range(0,8):
    sel=(omega[1:]==k); n=int(sel.sum())
    if n<200: om.append(dict(omega=k,n=n)); continue
    meas=float(tw[sel].mean()/p0); m2=float(mc[sel].mean()); mw=float(wm[sel].mean())
    om.append(dict(omega=k,n=n,ntwin=int(tw[sel].sum()),meas=meas,model_M2=m2,wind_M=mw))
    print("  w=%d n=%9d meas=%.4f  model<M2>=%.4f (rel %+.2f%%)  windM=%.4f"%(k,n,meas,m2,100*(meas/m2-1),mw))

# ---- also re-emit per-prime marginal (for the confirmed table) ----
pp=[]
for q in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]:
    d=np.zeros(NMAX+1,bool); d[q::q]=True; d[0]=False
    rate=is_twin[d].mean(); pp.append(dict(q=q,enrich=float(rate/p0),pred=q/(q-2)))

json.dump(dict(p0=p0,ntwin=ntwin,NMAX=NMAX,LIMIT=LIMIT,C=float(C),
   binned=binrows,per_omega=om,perprime=pp,
   primorial=dict(wind=wind,percentre=percentre,core=core)),
   open(os.path.join(DATA,"audit_percentre_results.json"),"w"),indent=1)
print("\nTOTAL %.1fs"%(time.time()-t0))
