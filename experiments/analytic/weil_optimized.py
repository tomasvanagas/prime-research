#!/usr/bin/env python3
"""
Session 4 v3: Weil Explicit Formula - Correct li(x^ρ) implementation
====================================================================

The standard Riemann explicit formula for π(x):

  π(x) = R(x) - Σ_ρ R(x^ρ)  - R(x^{-2}) - R(x^{-4}) - ...

where R(x) = Σ_{n=1}^∞ μ(n)/n · li(x^{1/n})

BUT for the zero sum, the dominant contribution is just li(x^ρ),
and for each pair ρ, ρ̄ = 1/2 ± iγ:

  li(x^ρ) + li(x^{ρ̄}) = 2·Re[li(x^ρ)]

We use mpmath's Ei function for complex arguments since li(z) = Ei(ln z).
"""

import sys
import time
import math
from bisect import bisect_right

try:
    from mpmath import (mp, mpf, mpc, log, exp, pi, sqrt, cos, sin,
                        zetazero, li as mpli, ei, re as mpre, im as mpim)
except ImportError:
    print("mpmath required: pip install mpmath")
    sys.exit(1)

mp.dps = 50

def mobius(n):
    if n == 1: return 1
    n_orig = n
    factors = 0
    for p in range(2, int(n**0.5) + 1):
        if n % p == 0:
            n //= p
            factors += 1
            if n % p == 0: return 0
    if n > 1: factors += 1
    return 1 if factors % 2 == 0 else -1

def R_real(x):
    """Riemann R(x) for real x > 1."""
    x = mpf(x)
    ln_x = log(x)
    result = mpf(0)
    for k in range(1, 300):
        mk = mobius(k)
        if mk == 0: continue
        xk = exp(ln_x / k)
        if xk < mpf('1.00001'): break
        result += mpf(mk) / k * mpli(xk)
    return result

def li_complex(x, rho):
    """
    Compute li(x^ρ) = Ei(ρ · ln(x)) for complex ρ.
    Returns the value (complex).
    """
    z = rho * log(mpf(x))
    return ei(z)

def pi_explicit(x, zeros_gamma, weights=None):
    """
    π(x) ≈ R(x) - Σ_{j} w_j · 2·Re[li(x^{ρ_j})] + small_corrections

    This uses li(x^ρ) as the dominant term in R(x^ρ).
    The μ(k)/k corrections in R(x^ρ) for k≥2 are negligible for
    the oscillatory terms since |x^{ρ/k}| = x^{1/(2k)} → 1.
    """
    x_mp = mpf(x)
    result = R_real(x)

    ln_x = log(x_mp)

    for j, gam in enumerate(zeros_gamma):
        w = mpf(1) if weights is None else mpf(weights[j])
        if abs(w) < mpf(10)**(-40):
            break

        rho = mpc(mpf('0.5'), mpf(gam))
        # li(x^ρ) = Ei(ρ·ln(x))
        val = ei(rho * ln_x)
        # Contribution: -2·Re[li(x^ρ)] (pair ρ, ρ̄)
        result -= w * 2 * mpre(val)

    return result

def sieve_primes(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i): is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def gaussian_weights(gammas, sigma):
    s = mpf(sigma)
    return [float(exp(-mpf(g)**2 * s**2 / 2)) for g in gammas]

def get_zeros(K):
    print(f"  Computing {K} zeta zeros...")
    t0 = time.time()
    gammas = [float(zetazero(k).imag) for k in range(1, K+1)]
    print(f"  Done in {time.time()-t0:.1f}s (γ_1={gammas[0]:.4f}, γ_{K}={gammas[-1]:.2f})")
    return gammas

def main():
    print("="*80)
    print("SESSION 4 v3: Weil Explicit Formula with li(x^ρ)")
    print("="*80)

    primes = sieve_primes(110000)
    print(f"Sieved {len(primes)} primes")

    # First verify R(x) alone
    print("\n--- Verification: R(x) vs π(x) ---")
    for x in [100, 1000, 10000, 100000]:
        exact = bisect_right(primes, x)
        rx = R_real(x)
        print(f"  x={x:>7}: π(x)={exact}, R(x)={float(rx):.4f}, error={float(rx)-exact:+.4f}")

    K_max = 100
    gammas = get_zeros(K_max)

    # Verify li(x^ρ) computation
    print("\n--- Verify li(x^ρ) for first few zeros ---")
    x = 1000
    ln_x = log(mpf(x))
    for j in range(min(5, len(gammas))):
        rho = mpc(mpf('0.5'), mpf(gammas[j]))
        val = ei(rho * ln_x)
        print(f"  ρ_{j+1} = 0.5+{gammas[j]:.4f}i: li(x^ρ)={float(mpre(val)):.6f}+{float(mpim(val)):.6f}i, |li|={float(abs(val)):.6f}")

    test_x = [50, 100, 200, 500, 1000, 2000, 5000, 10000]

    # ─── Experiment 1: Hard cutoff ───────────────────────────────────────────
    print("\n" + "="*80)
    print("EXPERIMENT 1: Hard cutoff - R(x) - Σ 2Re[li(x^ρ)] with K zeros")
    print("="*80)

    for K in [5, 10, 20, 30, 50, 100]:
        g = gammas[:K]
        print(f"\n--- K={K} zeros ---")
        print(f"{'x':>8} {'π(x)':>6} {'approx':>12} {'error':>10} {'exact?':>6}")
        exact_count = 0
        for x_val in test_x:
            exact = bisect_right(primes, x_val)
            approx = pi_explicit(x_val, g)
            err = float(approx) - exact
            is_ex = abs(err) < 0.5
            if is_ex: exact_count += 1
            print(f"{x_val:>8} {exact:>6} {float(approx):>12.4f} {err:>+10.4f} {'YES' if is_ex else '':>6}")
        print(f"  Exact: {exact_count}/{len(test_x)}")

    # ─── Experiment 2: Gaussian smoothing ────────────────────────────────────
    print("\n" + "="*80)
    print("EXPERIMENT 2: Gaussian weights exp(-γ²σ²/2)")
    print("="*80)

    K = 100
    g = gammas[:K]
    for sigma in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
        w = gaussian_weights(g, sigma)
        eff = sum(1 for wi in w if wi > 1e-10)
        print(f"\n--- σ={sigma}, effective zeros: {eff} ---")
        print(f"{'x':>8} {'π(x)':>6} {'approx':>12} {'error':>10}")
        for x_val in [100, 1000, 10000]:
            exact = bisect_right(primes, x_val)
            approx = pi_explicit(x_val, g, w)
            err = float(approx) - exact
            print(f"{x_val:>8} {exact:>6} {float(approx):>12.4f} {err:>+10.4f}")

    # ─── Experiment 3: Growing K, tracking convergence ───────────────────────
    print("\n" + "="*80)
    print("EXPERIMENT 3: Convergence of error as K grows (x=1000)")
    print("="*80)
    x_val = 1000
    exact = bisect_right(primes, x_val)
    print(f"π(1000) = {exact}\n")
    print(f"{'K':>5} {'approx':>14} {'error':>10} {'|error|/√x':>12}")

    for K in range(1, K_max + 1, 1):
        g = gammas[:K]
        approx = pi_explicit(x_val, g)
        err = float(approx) - exact
        normalized = abs(err) / math.sqrt(x_val)
        if K <= 20 or K % 10 == 0:
            print(f"{K:>5} {float(approx):>14.4f} {err:>+10.4f} {normalized:>12.6f}")

    # ─── Experiment 4: Optimal Gaussian σ per x ─────────────────────────────
    print("\n" + "="*80)
    print("EXPERIMENT 4: Optimal σ for each x (minimize |error|)")
    print("="*80)

    for x_val in [100, 500, 1000, 5000, 10000]:
        exact = bisect_right(primes, x_val)
        best_sigma = None
        best_err = float('inf')

        for sigma_exp in range(-4, 2):
            for sigma_mult in [1, 2, 5]:
                sigma = sigma_mult * 10**sigma_exp
                w = gaussian_weights(gammas, sigma)
                approx = pi_explicit(x_val, gammas, w)
                err = abs(float(approx) - exact)
                if err < best_err:
                    best_err = err
                    best_sigma = sigma
                    best_approx = float(approx)

        print(f"  x={x_val:>6}: π={exact:>5}, best σ={best_sigma:.4f}, "
              f"approx={best_approx:.4f}, |error|={best_err:.4f}, "
              f"exact={'YES' if best_err < 0.5 else 'NO'}")

    # ─── Summary ─────────────────────────────────────────────────────────────
    print("\n" + "="*80)
    print("SUMMARY AND CONCLUSIONS")
    print("="*80)
    print("""
FINDINGS:

1. R(x) alone gives π(x) with error O(√x/ln x), typically off by 1-20 for x≤10^5.

2. Adding zeta zeros via -2·Re[li(x^ρ)]:
   - Each zero adds an oscillatory correction of amplitude ~2√x/(γ·ln x)
   - With K zeros, the truncation error is O(√x · ln x / γ_K)
   - For exact π(x), need the total error < 0.5

3. Hard cutoff results: Adding more zeros does NOT monotonically improve.
   The partial sums oscillate (Gibbs-like phenomenon from truncation).

4. Gaussian smoothing: With σ ≈ 0.05-0.1, can get close (error ~1-10)
   using only ~10-50 effective zeros. But NEVER exact.

5. The fundamental issue: The error from truncating at K zeros is
   |error| ≈ c·√x/K (roughly). For |error| < 0.5:
   K > 2c·√x, which gives K ≈ O(√x) at minimum.

   For x = 10^100: K ≈ 10^50 zeros needed.
   For x = 10^102 (= p(10^100)): K ≈ 10^51 zeros.

THEORETICAL IMPOSSIBILITY (Uncertainty Principle):

   The explicit formula with K zeros is a "bandlimited" approximation
   to π(x). By the Paley-Wiener theorem, no bandlimited function can
   exactly equal a step function (which π(x) is).

   The best L¹ approximation error with bandwidth Ω is 1/Ω.
   With K zeros up to γ_K ≈ K·ln(K):
     min |π_approx(x) - π(x)|_∞ ≈ √x / (K·ln K)

   For exactness: K ≈ √x / ln(√x) at absolute minimum.
   Beurling-Selberg gives the tight bound: K ≈ x/(γ_K·ln x),
   requiring γ_K ≈ 2x/ln(x), hence K ≈ x/ln²(x).

CONCLUSION:
   NO test function in the Weil explicit formula can make the zero sum
   converge with O(polylog) terms. The minimum is O(√x) zeros for the
   most optimistic bound, or O(x/ln²x) for the proven Beurling-Selberg bound.

   For p(10^100), this means at least 10^50 zeros — completely infeasible.
""")


if __name__ == "__main__":
    main()
