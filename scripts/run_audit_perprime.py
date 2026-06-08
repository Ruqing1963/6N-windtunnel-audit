import numpy as np, time, json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); os.makedirs(DATA, exist_ok=True)
t0=time.time()
LIMIT = 200_000_000
NMAX  = LIMIT // 6           # = 100_000_000 centres

# ---- prime sieve to LIMIT+2 ----
sieve = np.ones(LIMIT+3, dtype=bool); sieve[:2]=False
for p in range(2, int((LIMIT+2)**0.5)+1):
    if sieve[p]:
        sieve[p*p::p]=False
print("prime sieve done %.1fs"%(time.time()-t0)); 

# ---- twin centre indicator for N=1..NMAX ----
# is_twin[N] = sieve[6N-1] & sieve[6N+1]
is_twin = np.zeros(NMAX+1, dtype=bool)
CH=10_000_000
for a in range(1, NMAX+1, CH):
    b=min(a+CH, NMAX+1)
    idx=np.arange(a,b,dtype=np.int64)
    is_twin[a:b] = sieve[6*idx-1] & sieve[6*idx+1]
ntwin=int(is_twin.sum())
print("twin centres up to 6N=%d : %d  (%.1fs)"%(LIMIT,ntwin,time.time()-t0))

# ---- omega_{>3}(N) for N=1..NMAX ----
omega=np.zeros(NMAX+1,dtype=np.int8)
primes_small = np.nonzero(sieve[:NMAX+1])[0]
primes_gt3 = primes_small[primes_small>=5]
for p in primes_gt3:
    omega[p::p]+=1
print("omega sieve done %.1fs  max omega=%d"%(time.time()-t0, int(omega.max())))

# baseline twin rate over centres N=1..NMAX
p0 = ntwin/NMAX
print("overall twin rate p0 = %.8f"%p0)

# ================= AUDIT 1: per-prime enrichment =================
# For prime q: rate among centres with q|N, and with q!|N.
res_perprime=[]
for q in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]:
    div = np.zeros(NMAX+1,dtype=bool); div[q::q]=True
    div[0]=False
    nd = div.sum()
    rate_div  = is_twin[div].mean()
    # complement over N=1..NMAX
    rate_ndiv = (ntwin - is_twin[div].sum())/(NMAX - nd)
    enrich_vs_all = rate_div/p0
    ratio_div_ndiv = rate_div/rate_ndiv
    res_perprime.append(dict(q=int(q),
        n_div=int(nd),
        rate_div=float(rate_div), rate_ndiv=float(rate_ndiv),
        enrich_vs_all=float(enrich_vs_all), pred_vs_all=q/(q-2),
        ratio_div_ndiv=float(ratio_div_ndiv), pred_div_ndiv=(q-1)/(q-3)))
    print("q=%2d  enrich(vs all)=%.5f pred q/(q-2)=%.5f | ratio(div/ndiv)=%.5f pred (q-1)/(q-3)=%.5f  (n_div=%d)"%(
        q, enrich_vs_all, q/(q-2), ratio_div_ndiv, (q-1)/(q-3), nd))

# ================= AUDIT 2: multiplicativity (pairwise & triple) =================
def rate_div_by(qs):
    m=np.zeros(NMAX+1,dtype=bool); m[1:]=True; m[0]=False
    for q in qs:
        d=np.zeros(NMAX+1,dtype=bool); d[q::q]=True
        m &= d
    n=m.sum()
    return is_twin[m].sum()/n, int(n)

res_mult=[]
combos=[(5,7),(5,11),(7,11),(5,13),(5,7,11),(5,7,13),(5,7,11,13)]
for qs in combos:
    r,n = rate_div_by(qs)
    meas = r/p0
    pred = 1.0
    for q in qs: pred*= q/(q-2)
    res_mult.append(dict(qs=list(map(int,qs)), n=n, meas=float(meas), pred=float(pred),
                         resid=float(meas/pred-1)))
    print("div by %-14s meas enrich=%.5f  prod pred=%.5f  resid=%+.3f%%  (n=%d)"%(
        str(qs), meas, pred, 100*(meas/pred-1), n))

# ================= AUDIT 3: omega-stratified enrichment =================
res_omega=[]
for k in range(0,8):
    sel = (omega[1:]==k)
    n=int(sel.sum())
    if n<50: 
        res_omega.append(dict(omega=k,n=n,enrich=None)); continue
    rate = is_twin[1:][sel].mean()
    res_omega.append(dict(omega=k,n=n,enrich=float(rate/p0)))
    print("omega=%d  n=%9d  enrich=%.4f"%(k,n,rate/p0))

out=dict(LIMIT=LIMIT,NMAX=NMAX,ntwin=ntwin,p0=p0,
         perprime=res_perprime,mult=res_mult,omega=res_omega,
         maxomega=int(omega.max()))
json.dump(out,open(os.path.join(DATA,"audit_perprime_results.json"),"w"),indent=1)
print("TOTAL %.1fs"%(time.time()-t0))
