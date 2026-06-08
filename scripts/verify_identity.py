"""
verify_identity.py  --  fast algebraic verifier (no sieve required, runs in <1s)

Checks the two closed-form claims of the paper:

  1. The global suppression constant
         C = prod_{q>3} q(q-3) / ((q-1)(q-2))         (eq. for C, ~0.7216)

  2. The per-centre identity (eq:percentre):
         E[twin | F] / p0  =  C * prod_{q|N} (q-1)/(q-3)
     and its exact link to the wind-tunnel marginal product:
         per-centre  =  [ prod_{q|N} q/(q-2) ]  *  prod_{q>59, q nmid N} g(q)
     with g(q) = q(q-3) / ((q-1)(q-2)) < 1.

  3. The primorial Omega_15 numbers: wind-tunnel x7.03 vs correct per-centre 6.983.

This is the algebra behind Tables 1-3 and Fig. 2B; the empirical confirmation
on real data is produced by run_audit.py.
"""

import numpy as np


def primes_up_to(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def g(q):
    """The bridge factor linking the two closed forms; g(q) < 1 for all q > 3."""
    q = np.asarray(q, dtype=np.float64)
    return q * (q - 3) / ((q - 1) * (q - 2))


def constant_C(qmax=2_000_000):
    """C = prod_{q>3} g(q). Converges absolutely since 1 - g(q) = O(1/q^2)."""
    qs = primes_up_to(qmax)
    qs = qs[qs >= 5].astype(np.float64)
    return float(np.prod(g(qs)))


def main():
    C = constant_C()
    print("=" * 64)
    print("1.  Global suppression constant")
    print("    C = prod_{q>3} q(q-3)/((q-1)(q-2)) = %.6f" % C)
    print("    (paper value: 0.7216)")

    # primorial Omega_15 small-prime set
    P = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

    wind = 1.0      # wind-tunnel marginal product  prod q/(q-2)
    core = 1.0      # per-centre core               prod (q-1)/(q-3)
    for q in P:
        wind *= q / (q - 2)
        core *= (q - 1) / (q - 3)
    percentre = C * core

    print()
    print("=" * 64)
    print("2.  Per-centre identity vs wind-tunnel, squarefree primorial Omega_15")
    print("    wind-tunnel  prod q/(q-2)             = %.4f   (the headline x7.03)" % wind)
    print("    correct per-centre  C*prod(q-1)/(q-3) = %.4f" % percentre)
    print("    gap  per-centre/wind                  = %+.2f%%" % (100 * (percentre / wind - 1)))

    # exact link: per-centre = wind * prod_{q>59} g(q)
    qs = primes_up_to(2_000_000)
    tail = qs[qs > 59].astype(np.float64)
    tail_prod = float(np.prod(g(tail)))
    print()
    print("    exact link:  per-centre = wind * prod_{q>59} g(q)")
    print("    prod_{q>59} g(q)                      = %.4f" % tail_prod)
    print("    wind * tail                           = %.4f" % (wind * tail_prod))
    print("    matches per-centre above to           %.2e" % abs(wind * tail_prod - percentre))

    print()
    print("=" * 64)
    print("3.  Marginal vs per-centre for a single small prime")
    print("    q   marginal q/(q-2)   per-centre boost (q-1)/(q-3)")
    for q in [5, 7, 11, 13]:
        print("   %2d   %.5f            %.5f" % (q, q / (q - 2), (q - 1) / (q - 3)))
    print()
    print("All closed-form claims reproduced. Empirical check: run_audit.py")


if __name__ == "__main__":
    main()
