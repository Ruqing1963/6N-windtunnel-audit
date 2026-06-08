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
# omega and model multiplier M(N)=prod_{q|N,q>3} q/(q-2)
omega=np.zeros(NMAX+1,dtype=np.int8)
M=np.ones(NMAX+1,dtype=np.float64)
primes=np.nonzero(sieve[:NMAX+1])[0]; primes=primes[primes>=5]
for q in primes:
    omega[q::q]+=1; M[q::q]*= q/(q-2)
print("built omega,M  %.1fs"%(time.time()-t0))

# per-omega: measured enrichment vs model <M> vs primorial envelope
def prim_env(k):
    P=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]; v=1.0
    for i in range(k): v*=P[i]/(P[i]-2)
    return v
rows=[]
for k in range(0,9):
    sel=(omega[1:]==k); n=int(sel.sum())
    if n<30: rows.append(dict(omega=k,n=n)); continue
    meas=float(is_twin[1:][sel].mean()/p0)
    modM=float(M[1:][sel].mean())          # correct multiplicative prediction (avg over real factorisations)
    ntw=int(is_twin[1:][sel].sum())
    rows.append(dict(omega=k,n=n,ntwin=ntw,meas=meas,model_M=modM,
                     primorial=prim_env(k),rel_err_modelM=meas/modM-1))
    print("w=%d n=%9d ntw=%8d  meas=%.4f  <M>model=%.4f (rel %+.2f%%)  primorial=%.4f"%(
        k,n,ntw,meas,modM,100*(meas/modM-1),prim_env(k)))

# global multiplicative-model check: E[is_twin]/p0 vs M, correlation
# slope of is_twin/p0 on M should be ~1 if model exact
mean_M_all=float(M[1:].mean()); mean_M_twin=float(M[1:][is_twin[1:]].mean())
print("mean M over all=%.5f  over twins=%.5f  ratio=%.5f (pred ~ E[M^2]/E[M] enrichment)"%(
    mean_M_all,mean_M_twin,mean_M_twin/mean_M_all))

# binned by model multiplier M: does measured enrichment track M 1:1?
edges=np.array([0.9,1.05,1.2,1.4,1.6,1.9,2.3,3.0,5.0])
Mc=M[1:]; tw=is_twin[1:]
binrows=[]
for i in range(len(edges)-1):
    sel=(Mc>=edges[i])&(Mc<edges[i+1]); n=int(sel.sum())
    if n<1000: continue
    meas=float(tw[sel].mean()/p0); mid=float(Mc[sel].mean())
    binrows.append(dict(lo=float(edges[i]),hi=float(edges[i+1]),n=n,meanM=mid,meas=meas,ratio=meas/mid))
    print("M in [%.2f,%.2f) n=%9d  <M>=%.4f  meas_enrich=%.4f  ratio=%.4f"%(edges[i],edges[i+1],n,mid,meas,meas/mid))

# primorial envelope out to omega=15 (the wind-tunnel curve) -> x7.03 at 15
env=[prim_env(k) for k in range(0,16)]
print("primorial envelope omega=0..15:", [round(e,3) for e in env], "  -> x%.3f at 15"%env[15])

json.dump(dict(p0=p0,ntwin=ntwin,NMAX=NMAX,LIMIT=LIMIT,
   per_omega=rows,mean_M_all=mean_M_all,mean_M_twin=mean_M_twin,
   binned=binrows,prim_env=env),open(os.path.join(DATA,"audit_omega_results.json"),"w"),indent=1)
print("TOTAL %.1fs"%(time.time()-t0))
