# 6N Wind-Tunnel Audit

**Paper 28 / Part XXVIII** of the *Arithmetic Geodynamics on the 6N Skeleton* series.

Auditing the Part XXIII "wind-tunnel" extrapolation of twin-centre enrichment
($\times 7.03$ at the $\omega=15$ primorial) against a **real prime sieve** to
$6N \le 2\times10^{8}$ (813,370 twin centres). No fits, no fabricated data — every
number in the paper is reproduced by the scripts here.

> A *twin centre* is an integer $N$ for which both $6N-1$ and $6N+1$ are prime.
> $\omega_{>3}(N)$ counts the distinct prime factors of $N$ that exceed 3.

---

## The verdict in one paragraph

The wind-tunnel **mechanism is vindicated**: the per-prime modifier $q/(q-2)$ is
exact to $\lesssim 0.3\%$ for every $q = 5,\dots,59$, and its multiplicativity holds
within Poisson error across the reachable range. The headline $\times 7.03$ is
**sound for what it actually describes** — the squarefree primorial centre (to 0.66%)
and the marginal divisibility statement. The audit adds two clarifications and no new
theorem:

1. **Marginal vs. per-centre (interpretive).** $\prod_{q\mid N} q/(q-2)$ is a *marginal*
   factor. The correct *per-centre* likelihood is
   $$\mathbb{E}[\text{twin}\mid N]/p_0 \;=\; \mathcal{C}\prod_{q\mid N}\frac{q-1}{q-3},
   \qquad \mathcal{C}=\prod_{q>3}\frac{q(q-3)}{(q-1)(q-2)}=0.7216,$$
   which tracks the data to $<1\%$ everywhere. Using the bare wind-tunnel product as a
   per-centre multiplier over-predicts non-primorial factorisations by up to ~38%.
2. **Envelope vs. typical (scope).** $\times 7.03$ is an *upper envelope* realised only by
   the vanishingly rare primorial-loaded centre. The typical $\omega$-stratified average
   enrichment sits far below it (measured $\langle E\rangle\approx 2.63$ at $\omega=6$ vs.
   the 4.27 envelope).

This is a measurement, not a theorem; no claim is made about the infinitude of any
constellation.

---

## Repository layout

```
.
├── paper/
│   ├── paper28.pdf                     # the compiled paper (6 pp.)
│   └── paper28.tex                     # LaTeX source
├── code/
│   ├── verify_identity.py              # fast (<1s) algebraic verifier — no sieve needed
│   ├── run_audit_perprime.py           # per-prime law + multiplicativity (Tables 1, 2)
│   ├── run_audit_omega.py              # omega-stratified + M-binned (Table 3, Fig 2A)
│   ├── run_audit_percentre.py          # per-centre model + primorial (Fig 1B, Fig 2B)
│   └── make_figures.py                 # regenerates both figure PDFs from the data
├── data/
│   ├── audit_perprime_results.json
│   ├── audit_omega_results.json
│   └── audit_percentre_results.json
├── figures/
│   ├── fig1_audit.pdf
│   └── fig2_envelope.pdf
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

---

## Reproducing the results

Requirements: Python 3.9+, `numpy`, `matplotlib` (see `requirements.txt`).

```bash
pip install -r requirements.txt
```

**Fast check (no sieve, < 1 second)** — reproduces the constant $\mathcal{C}=0.7216$, the
$7.03$-vs-$6.983$ primorial gap, and the exact link between the two closed forms:

```bash
python code/verify_identity.py
```

**Full empirical audit (~35 s total, ~2 GB RAM).** Each script is self-contained: it
sieves all primes to $2\times10^{8}$, marks the 813,370 twin centres, and writes its JSON
into `data/`. Run them in any order, then build the figures:

```bash
python code/run_audit_perprime.py     # ~7 s
python code/run_audit_omega.py        # ~12 s
python code/run_audit_percentre.py    # ~16 s
python code/make_figures.py           # regenerates figures/*.pdf from data/
```

The regenerated JSON is bit-for-bit identical to the committed files, and the figures
match those in the paper.

### Hardware note
The sieve allocates a boolean array of length $2\times10^{8}$ plus per-centre `omega`/`M`
arrays over $N\le 3.3\times10^{7}$; peak usage is roughly 2 GB. To run on a smaller
machine, lower `LIMIT` at the top of each script (e.g. `LIMIT=60_000_000`); the closed-form
agreements hold at any scale, only the Poisson error bars widen.

---

## What each output contains

- **`audit_perprime_results.json`** — for each $q\in\{5,\dots,59\}$: measured marginal
  enrichment vs. $q/(q-2)$, divisor/non-divisor ratio vs. $(q-1)/(q-3)$; pairwise/triple
  multiplicativity residuals; raw $\omega$-stratified rates.
- **`audit_omega_results.json`** — $\omega$-stratified measured enrichment vs. the
  average multiplicative model $\langle M\rangle$ and the primorial envelope; $M$-binned
  ratios; the primorial envelope curve out to $\omega=15$.
- **`audit_percentre_results.json`** — the constant $\mathcal{C}$; per-centre model
  $\mathcal{C}\prod(q-1)/(q-3)$ binned against data; per-$\omega$ comparison; the primorial
  $\Omega_{15}$ numbers (`wind` $=7.0299$, `percentre` $=6.9833$, `core`).

---

## The series

This is Part XXVIII. Directly relevant predecessors:

- **Part XXIII** — *Phase transition kinetics and topological event horizons at extreme
  $\omega$*, [doi:10.5281/zenodo.20586919](https://doi.org/10.5281/zenodo.20586919)
  (the source of the $\times 7.03$ figure audited here).
- **Part XXVI** — *The second moment of $\omega$ under twin-centre conditioning on the 6N
  skeleton*, [doi:10.5281/zenodo.20593654](https://doi.org/10.5281/zenodo.20593654)
- **Review of Parts I–XIX**,
  [doi:10.5281/zenodo.20585301](https://doi.org/10.5281/zenodo.20585301)

---

## Citation

See `CITATION.cff`, or:

> R. Chen, *Auditing the Wind Tunnel: the Twin-Enrichment Closed Form at Reachable
> $\omega$, and the Marginal-versus-Per-Centre Distinction* (Part XXVIII), 2026.
> https://github.com/Ruqing1963/6N-windtunnel-audit

## License

Code under the MIT License (`LICENSE`). The paper text and figures are © 2026 Ruqing Chen.
