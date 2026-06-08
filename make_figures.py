import json, numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,
    "figure.dpi":150,"savefig.bbox":"tight"})

A=json.load(open(os.path.join(DATA,"audit_perprime_results.json")))   # per-prime (both forms), mult
A3=json.load(open(os.path.join(DATA,"audit_percentre_results.json")))  # binned M2/Mwind, per_omega model, primorial, C
A2=json.load(open(os.path.join(DATA,"audit_omega_results.json")))  # prim_env

C=A3["C"]; p0=A3["p0"]

# ---------- FIGURE 1 : the audit ----------
fig,ax=plt.subplots(1,2,figsize=(9.2,4.0))

# (A) per-prime marginal law confirmed
pp=A["perprime"]
q=np.array([r["q"] for r in pp])
e_meas=np.array([r["enrich_vs_all"] for r in pp]); e_pred=np.array([r["pred_vs_all"] for r in pp])
r_meas=np.array([r["ratio_div_ndiv"] for r in pp]); r_pred=np.array([r["pred_div_ndiv"] for r in pp])
ax[0].plot([1,1.75],[1,1.75],'k-',lw=0.8,alpha=0.6,label="exact (y=x)")
ax[0].scatter(e_pred,e_meas,s=28,c="#1f77b4",zorder=3,label=r"marginal $q/(q-2)$")
ax[0].scatter(r_pred,r_meas,s=28,marker="s",facecolors='none',edgecolors="#d62728",zorder=3,
              label=r"div/non-div $(q-1)/(q-3)$")
for qi,xp,yp in zip(q,e_pred,e_meas):
    if qi in (5,7,11): ax[0].annotate("q=%d"%qi,(xp,yp),textcoords="offset points",xytext=(4,-9),fontsize=7)
ax[0].set_xlabel("closed-form prediction"); ax[0].set_ylabel("measured enrichment")
ax[0].set_title("(A) microscopic law confirmed (q=5..59)",fontsize=9)
ax[0].legend(fontsize=7,loc="upper left")
# residual inset
axin=ax[0].inset_axes([0.56,0.12,0.40,0.34])
axin.axhline(0,color='k',lw=0.6)
axin.scatter(q,100*(e_meas/e_pred-1),s=10,c="#1f77b4")
axin.set_title(r"$100(\mathrm{meas}/\mathrm{pred}-1)$",fontsize=6); axin.set_xlabel("q",fontsize=6)
axin.tick_params(labelsize=6); axin.set_ylim(-1,1)

# (B) per-centre reconciliation
b=A3["binned"]
M2=np.array([r["M2"] for r in b]); Mw=np.array([r["Mwind"] for r in b]); meas=np.array([r["meas"] for r in b])
lo=min(M2.min(),Mw.min(),meas.min())-0.1; hi=max(Mw.max(),meas.max())+0.2
ax[1].plot([lo,hi],[lo,hi],'k-',lw=0.8,alpha=0.6,label="exact (y=x)")
ax[1].scatter(M2,meas,s=30,c="#2ca02c",zorder=3,
   label=r"vs correct  $C\!\prod\frac{q-1}{q-3}$")
ax[1].scatter(Mw,meas,s=30,marker="^",c="#ff7f0e",zorder=3,
   label=r"vs wind-tunnel $\prod\frac{q}{q-2}$")
# primorial markers
pw=A3["primorial"]["wind"]; pc=A3["primorial"]["percentre"]
ax[1].set_xlabel("model per-centre enrichment"); ax[1].set_ylabel("measured per-centre enrichment")
ax[1].set_title(r"(B) per-centre: wind-tunnel $M$ overshoots; $C\!=\!%.4f$"%C,fontsize=9)
ax[1].legend(fontsize=7,loc="upper left")
fig.suptitle("Auditing the twin-enrichment closed form against real centres (6N $\\leq$ 2$\\times$10$^8$, %s twin centres)"%
             format(A3["ntwin"],","),fontsize=9.5)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(os.path.join(FIG,"fig1_audit.pdf")); print("fig1 done")

# ---------- FIGURE 2 : envelope correction ----------
fig2,ax2=plt.subplots(1,2,figsize=(9.2,3.9))
om=[r for r in A3["per_omega"] if r.get("meas") is not None and r["omega"]>=1]
ks=np.array([r["omega"] for r in om]); ms=np.array([r["meas"] for r in om])
mod=np.array([r["model_M2"] for r in om]); ns=np.array([r["n"] for r in om]); nt=np.array([r["ntwin"] for r in om])
env=A2["prim_env"]  # index = omega 0..15
kk=np.arange(0,16)
ax2[0].plot(kk,env,'-',color="#9467bd",lw=1.6,label=r"primorial envelope $\prod_{q\leq p_\omega}\frac{q}{q-2}$")
ax2[0].scatter([15],[env[15]],c="#9467bd",s=40,zorder=4)
ax2[0].annotate(r"$\times7.03$ at $\omega=15$ (Part XXIII)",(15,env[15]),
                textcoords="offset points",xytext=(-150,-2),fontsize=7,color="#9467bd")
yerr=ms/np.sqrt(np.maximum(nt,1))
ax2[0].errorbar(ks,ms,yerr=yerr,fmt='o',ms=5,c="#1f77b4",capsize=2,zorder=3,
                label=r"measured $\langle$enrich$\rangle(\omega)$ (data)")
ax2[0].plot(ks,mod,'s--',ms=4,c="#2ca02c",label=r"correct model $\langle C\!\prod\frac{q-1}{q-3}\rangle$")
ax2[0].set_xlabel(r"$\omega_{>3}(N)$"); ax2[0].set_ylabel("twin enrichment vs all centres")
ax2[0].set_title("(A) the wind-tunnel tilt: data sits far below the\nprimorial envelope at matched $\\omega$",fontsize=8.5)
ax2[0].legend(fontsize=7,loc="upper left"); ax2[0].set_xlim(0.5,15.5)

# (B) primorial centre: wind vs correct nearly coincide
labels=["wind-tunnel\n$\\prod q/(q-2)$","correct per-centre\n$C\\prod(q-1)/(q-3)$"]
vals=[pw,pc]
bars=ax2[1].bar([0,1],vals,color=["#ff7f0e","#2ca02c"],width=0.55)
ax2[1].set_xticks([0,1]); ax2[1].set_xticklabels(labels,fontsize=7.5)
ax2[1].set_ylabel("enrichment of the squarefree primorial $\\Omega_{15}$")
for x,v in zip([0,1],vals): ax2[1].text(x,v+0.06,"%.3f"%v,ha="center",fontsize=9)
ax2[1].set_ylim(0,7.8)
ax2[1].set_title("(B) for the primorial centre the two agree\nto %.2f%% (7.03 is sound there)"%(100*(pc/pw-1)),fontsize=8.5)
fig2.suptitle("The $\\times7.03$ figure: correct as a marginal / primorial number, but an upper envelope vs typical $\\omega$",fontsize=9.5)
fig2.tight_layout(rect=[0,0,1,0.95])
fig2.savefig(os.path.join(FIG,"fig2_envelope.pdf")); print("fig2 done")
