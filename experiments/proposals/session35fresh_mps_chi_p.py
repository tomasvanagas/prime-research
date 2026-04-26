"""
Proposal 3 — MPS / tensor-network bond dimension of chi_P along the natural
midpoint cut.

For L = 6..14, build chi_P[0..2^L), reshape into 2^{L/2} x 2^{L/2} matrix,
SVD, and report:
  - singular value decay (raw + ratios)
  - effective rank at tolerances 1e-2, 1e-4, 1e-6
  - Renyi-2 entropy
  - rank growth rate vs L

VERDICT: effective-rank-at-1e-4 stays below ~64 as L increases ⇒ escalate.
"""
import math
import numpy as np
from sympy import isprime


def build_chi(N):
    chi = np.zeros(N, dtype=np.float64)
    for k in range(2, N):
        if isprime(k):
            chi[k] = 1.0
    return chi


def schmidt_singular_values(chi):
    n = len(chi)
    L = int(round(math.log2(n)))
    half = L // 2
    rows = 1 << half
    cols = n // rows
    M = chi.reshape(rows, cols)
    sv = np.linalg.svd(M, compute_uv=False)
    return sv


def effective_rank(sv, tol):
    if sv[0] <= 0:
        return 0
    norm = sv / sv[0]
    return int(np.sum(norm > tol))


def renyi2_entropy(sv):
    s2 = sv ** 2
    p = s2 / s2.sum()
    p = p[p > 0]
    H2 = -math.log2(float((p ** 2).sum()))
    return H2


def main():
    print("# MPS / Schmidt-decomposition along midpoint cut of chi_P")
    print(f"# Each row: L, N=2^L, midpoint cut, top-5 singular values, eff_rank @1e-2/1e-4/1e-6, Renyi2")

    eff_at_1e4 = []
    Ls = list(range(6, 15))
    for L in Ls:
        N = 1 << L
        chi = build_chi(N)
        sv = schmidt_singular_values(chi)
        top5 = ", ".join(f"{s:.4f}" for s in sv[:5])
        e2 = effective_rank(sv, 1e-2)
        e4 = effective_rank(sv, 1e-4)
        e6 = effective_rank(sv, 1e-6)
        H2 = renyi2_entropy(sv)
        eff_at_1e4.append(e4)
        print(f"L={L:>2}  N={N:>5}  sv[0..4]=[{top5}]  "
              f"eff_rank(1e-2)={e2:>3}  (1e-4)={e4:>3}  (1e-6)={e6:>3}  "
              f"H2={H2:.3f}")

    # Fit eff_rank(1e-4) vs L
    arr = np.array(eff_at_1e4, dtype=float)
    Ls_arr = np.array(Ls, dtype=float)
    # Try: eff_rank ~ a * L^b vs eff_rank ~ a * 2^(b*L)
    # Log-log fit (against L)
    if (arr > 0).all():
        loglog = np.polyfit(np.log(Ls_arr), np.log(arr), 1)
        print(f"\n## eff_rank(1e-4) ~ L^{loglog[0]:.3f}  (intercept exp={math.exp(loglog[1]):.3f})")
        # Linear fit log(arr) ~ alpha * L (exponential growth check)
        lin = np.polyfit(Ls_arr, np.log2(arr), 1)
        print(f"## eff_rank(1e-4) ~ 2^({lin[0]:.3f} L)  ⇒ "
              f"polylog if alpha << 0.5; volume-law if alpha ≈ 0.5")
        if lin[0] < 0.2:
            print("\nVERDICT: eff-rank growth is polylog-compatible ⇒ "
                  "PROPOSAL 3 PROMISING (escalate to permuted orderings)")
        elif lin[0] < 0.4:
            print("\nVERDICT: sub-volume-law but worse than polylog "
                  "⇒ AMBIGUOUS")
        else:
            print("\nVERDICT: near-volume-law midpoint entanglement "
                  "⇒ PROPOSAL 3 FAILS at natural cut")

    print("\n# DONE")


if __name__ == "__main__":
    main()
