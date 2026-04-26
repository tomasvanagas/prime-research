"""
Walsh-Hadamard L1 Spectral Norm of the Prime Indicator
=======================================================

Hypothesis: If the prime indicator chi_P : {0,1}^k -> {0,1} has small
spectral L1 norm

    || chi_P_hat ||_1  =  sum_S | chi_P_hat(S) |

then by the Kushilevitz-Mansour algorithm, chi_P is computable by a
polynomial-size decision tree (in particular has a polylog approximate
representation), and pi(x) on x = 2^k can be evaluated in time
poly(L1(chi_P), k, 1/eps). A polylog L1 norm => polylog pi(x).

Mansour's theorem: any function with L1 spectral norm <= M is
eps-approximable by an O(M^2/eps^2)-sparse polynomial in the Walsh basis.

Prior project work measured Fourier WEIGHT by degree (W_k =
sum_{|S|=k} chi_P_hat(S)^2), which is degree-distribution but does NOT
upper-bound L1. The key inequality is

    || f_hat ||_1  >=  || f_hat ||_2  =  || f ||_2

so for chi_P with density rho, L1 >= sqrt(rho) >= 1/sqrt(log N).
Random functions of density rho have L1 ~ 2^{k/2} * sqrt(rho) (Khintchine).
A polylog L1 would be a HUGE structural anomaly.

This experiment computes:
  * L1 = sum |chi_P_hat(S)|         (Mansour-relevant)
  * L0 = #{ S : |chi_P_hat(S)| > tau }  (sparsity)
  * Linfty = max |chi_P_hat(S)|     (max correlation)

across k = 8..18, comparing to a random control of matching density.

Reads: chi_P viewed as boolean function on {0, 1, ..., 2^k - 1} indexed
by the natural binary encoding.
"""

import numpy as np
import time
from sympy import sieve


def prime_indicator_vec(k):
    N = 1 << k
    primes = list(sieve.primerange(2, N))
    v = np.zeros(N, dtype=np.float64)
    for p in primes:
        v[p] = 1.0
    return v


def random_indicator_vec(k, density, seed=0):
    N = 1 << k
    rng = np.random.default_rng(seed)
    n_one = max(1, int(round(density * N)))
    idx = rng.choice(N, size=n_one, replace=False)
    v = np.zeros(N, dtype=np.float64)
    v[idx] = 1.0
    return v


def walsh_hadamard_inplace(v):
    """In-place WHT (Sylvester ordering). Returns hat values WITHOUT 1/N normalization."""
    a = v.copy()
    n = len(a)
    h = 1
    while h < n:
        # Vectorized butterfly
        for i in range(0, n, h * 2):
            x = a[i:i + h].copy()
            y = a[i + h:i + 2 * h].copy()
            a[i:i + h] = x + y
            a[i + h:i + 2 * h] = x - y
        h *= 2
    return a / n  # convention: f_hat(S) = (1/N) sum f(x) (-1)^{<S,x>}


def measure_walsh(v):
    hat = walsh_hadamard_inplace(v)
    abshat = np.abs(hat)
    L1 = float(abshat.sum())
    L2 = float(np.sqrt((abshat ** 2).sum()))
    Linf = float(abshat.max())
    # L0 at fractional thresholds of Linf
    thresh = [Linf / 2, Linf / 10, Linf / 100, Linf / 1000]
    L0_at = {t: int((abshat > t).sum()) for t in thresh}
    # Also: total # nonzero (above 1e-15 numeric noise)
    L0_nonzero = int((abshat > 1e-15).sum())
    return {
        "L1": L1,
        "L2": L2,
        "Linf": Linf,
        "L1_over_L2": L1 / L2 if L2 > 0 else float("nan"),
        "L0_thresholds": L0_at,
        "L0_nonzero": L0_nonzero,
    }


def main():
    print("=" * 88)
    print("Walsh-Hadamard L1 Spectral Norm of chi_P")
    print("=" * 88)
    print()
    print(f"{'k':>3} {'N':>10} {'rho':>8} {'L1(P)':>10} {'L1(R)':>10} "
          f"{'L1(P)/sqrt(N)':>14} {'L1(R)/sqrt(N)':>14} {'L1(P)/L2(P)':>12}")
    print("-" * 88)

    rows = []
    for k in range(8, 19):
        N = 1 << k
        v_p = prime_indicator_vec(k)
        rho = float(v_p.sum() / N)
        v_r = random_indicator_vec(k, rho, seed=1234 + k)
        m_p = measure_walsh(v_p)
        m_r = measure_walsh(v_r)

        # Predicted scaling for random of density rho:
        # E[ |hat(S)| ] ~ sqrt(rho/N), so L1 ~ N * sqrt(rho/N) = sqrt(rho * N)
        # That's sqrt(rho) * sqrt(N).
        pred_random_L1_over_sqrtN = np.sqrt(rho)
        rows.append({
            "k": k, "N": N, "rho": rho,
            "L1_P": m_p["L1"], "L1_R": m_r["L1"],
            "Linf_P": m_p["Linf"], "Linf_R": m_r["Linf"],
            "L1_over_L2_P": m_p["L1_over_L2"],
            "L0_nonzero_P": m_p["L0_nonzero"],
        })
        print(f"{k:>3} {N:>10} {rho:>8.4f} "
              f"{m_p['L1']:>10.4f} {m_r['L1']:>10.4f} "
              f"{m_p['L1'] / np.sqrt(N):>14.4f} {m_r['L1'] / np.sqrt(N):>14.4f} "
              f"{m_p['L1_over_L2']:>12.4f}")

    print()
    print("Predicted random scaling: L1 / sqrt(N) ~ sqrt(rho) ~ sqrt(1/log N)")
    print()

    # --- Scaling fit: L1(P) = C * N^alpha
    print("Power-law fit: L1 = C * N^alpha (over k>=10)")
    sub = [r for r in rows if r["k"] >= 10]
    logN = np.log(np.array([r["N"] for r in sub]))
    logL1_P = np.log(np.array([r["L1_P"] for r in sub]))
    logL1_R = np.log(np.array([r["L1_R"] for r in sub]))
    A = np.vstack([logN, np.ones_like(logN)]).T
    a_p, c_p = np.linalg.lstsq(A, logL1_P, rcond=None)[0]
    a_r, c_r = np.linalg.lstsq(A, logL1_R, rcond=None)[0]
    print(f"  primes: alpha = {a_p:.4f}  (polylog -> 0; volume -> 0.5)")
    print(f"  random: alpha = {a_r:.4f}  (expected ~0.5 modulo log(rho))")

    # --- L0 sparsity at largest k
    print()
    last_k = max(r["k"] for r in rows)
    v_p = prime_indicator_vec(last_k)
    rho = float(v_p.sum() / (1 << last_k))
    v_r = random_indicator_vec(last_k, rho, seed=999)
    m_p = measure_walsh(v_p)
    m_r = measure_walsh(v_r)
    print(f"At k = {last_k} (N = {1 << last_k}):")
    print(f"  Linf(P) = {m_p['Linf']:.6e}   Linf(R) = {m_r['Linf']:.6e}")
    print(f"  Spectral L0 sparsity at fractional thresholds of Linf:")
    print(f"    threshold     |  primes count  |  random count")
    for t_lab, t in [("Linf/2", m_p['Linf']/2), ("Linf/10", m_p['Linf']/10),
                     ("Linf/100", m_p['Linf']/100), ("Linf/1000", m_p['Linf']/1000)]:
        # recompute on both with the SAME absolute threshold for comparability
        hatP = np.abs(walsh_hadamard_inplace(v_p))
        hatR = np.abs(walsh_hadamard_inplace(v_r))
        nP = int((hatP > t).sum())
        nR = int((hatR > t).sum())
        print(f"    {t_lab:10}     | {nP:>14} | {nR:>13}")

    # --- Verdict
    print()
    print("Interpretation:")
    if a_p < 0.1:
        print("  L1(P) sub-polynomial in N. INVESTIGATE -- could yield polylog algorithm.")
    elif abs(a_p - a_r) < 0.05:
        print("  L1(P) scales same as random. No exploitable structure for "
              "Kushilevitz-Mansour-style fast computation.")
    else:
        print(f"  L1(P) - L1(R) exponent gap = {a_r - a_p:.4f}: marginal.")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.2f}s")
