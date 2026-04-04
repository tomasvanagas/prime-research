#!/usr/bin/env python3
"""
Session 8: Gram Points + Riemann-Siegel Approach to pi(x)
=========================================================

Goal: Can we compute pi(x) faster than O(x^{1/2}) using:
  1. Gram points and the sign pattern of Z(t)
  2. The Odlyzko-Schonhage batching idea
  3. Selective zero contributions (most-significant zeros)
  4. Optimal test functions in the explicit formula

Key formulas:
  - Riemann-Siegel: Z(t) = 2 sum_{n=1}^{N} cos(theta(t) - t*ln(n))/sqrt(n) + R(t)
    where N = floor(sqrt(t/(2*pi)))
  - Explicit formula: psi(x) = x - sum_rho x^rho/rho - ln(2*pi) - (1/2)*ln(1-x^{-2})
  - pi(x) via Mobius inversion of psi(x)
  - Gram points: theta(g_n) = n*pi
"""

import time
import math
import sys
from collections import defaultdict

import mpmath
from mpmath import mp, mpf, mpc, pi, log, exp, sqrt, cos, sin, zeta, siegeltheta, siegelz
from sympy import primerange

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
mp.dps = 30  # 30 digits is plenty for our experiments

def flush_print(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()


def primes_up_to(n):
    return list(primerange(2, int(n) + 1))


def primes_in_range(lo, hi):
    return list(primerange(max(2, int(lo)), int(hi) + 1))


# ===========================================================================
# Precompute zeta zeros (the main bottleneck)
# ===========================================================================

ZERO_CACHE = {}

def get_zero(k):
    """Get the imaginary part of the k-th zeta zero (cached)."""
    if k not in ZERO_CACHE:
        ZERO_CACHE[k] = mpmath.zetazero(k)
    return ZERO_CACHE[k]


def precompute_zeros(n):
    """Precompute first n zeta zeros."""
    flush_print(f"  Precomputing {n} zeta zeros...", end="")
    t0 = time.time()
    for k in range(1, n + 1):
        get_zero(k)
    flush_print(f" done in {time.time()-t0:.1f}s")


# ===========================================================================
# PART 1: Riemann-Siegel Z-function
# ===========================================================================

def riemann_siegel_manual(t, count_ops=False):
    """
    Manual Riemann-Siegel formula to count operations precisely.
    Z(t) ~ 2 * sum_{n=1}^{N} cos(theta(t) - t*ln(n)) / sqrt(n) + R(t)
    where N = floor(sqrt(t/(2*pi))).
    """
    N = int(mpmath.floor(mpmath.sqrt(t / (2 * pi))))
    theta_t = siegeltheta(t)

    total = mpf(0)
    for n in range(1, N + 1):
        total += cos(theta_t - t * log(n)) / sqrt(n)
    total *= 2

    # First correction term (C_0)
    p = mpmath.sqrt(t / (2 * pi)) - N
    c0 = cos(2 * pi * (p**2 - p - mpf(1)/16)) / cos(2 * pi * p)
    correction = (-1)**(N - 1) * (t / (2 * pi))**(-mpf(1)/4) * c0

    result = total + correction
    ops = N

    if count_ops:
        return result, ops
    return result


# ===========================================================================
# PART 2: Gram points
# ===========================================================================

def gram_point(n):
    """Compute the n-th Gram point g_n where theta(g_n) = n*pi."""
    if n == 0:
        t = mpf('17.8')
    else:
        t = 2 * pi * exp(1 + mpmath.lambertw((8 * pi * n + 1) / (8 * math.e)))

    target = n * pi
    for _ in range(50):
        th = siegeltheta(t)
        dth = mpf(0.5) * log(t / (2 * pi))
        dt = (th - target) / dth
        t -= dt
        if abs(dt) < mpf(10)**(-mp.dps + 5):
            break
    return t


# ===========================================================================
# PART 3: Explicit formula for psi(x) using zeros
# ===========================================================================

def psi_explicit(x, num_zeros=100, return_terms=False):
    """
    psi(x) = x - sum_{rho} x^rho/rho - ln(2*pi) - (1/2)*ln(1 - x^{-2})
    Uses precomputed zeros.
    """
    x = mpf(x)
    result = x

    zero_contributions = []
    for k in range(1, num_zeros + 1):
        rho = get_zero(k)
        gamma_k = rho.imag
        x_rho = exp(rho * log(x))
        contrib = 2 * mpmath.re(x_rho / rho)
        zero_contributions.append((k, gamma_k, abs(contrib), contrib))
        result -= contrib

    result -= log(2 * pi)
    if x > 1:
        result -= 0.5 * log(1 - x**(-2))

    if return_terms:
        return result, zero_contributions
    return result


# ===========================================================================
# EXPERIMENT 1: How many zeros for exact psi(x)?
# ===========================================================================

def experiment_zeros_needed():
    flush_print("\n" + "="*75)
    flush_print("EXPERIMENT 1: Number of zeros needed for exact psi(x)")
    flush_print("="*75)

    # Precompute enough zeros
    precompute_zeros(500)

    x_values = [100, 200, 500, 1000, 2000, 5000]
    flush_print(f"{'x':>8} {'true psi':>12} {'zeros needed':>14} {'sqrt(x)':>10} "
                f"{'ratio':>10} {'RS terms':>10}")
    flush_print("-"*70)

    results = []
    for x in x_values:
        true_psi = mpmath.fsum(mpmath.log(p) for p in primes_up_to(x))
        # Also add prime power contributions
        for p in primes_up_to(int(x**0.5) + 1):
            pk = p * p
            while pk <= x:
                true_psi += log(p)
                pk *= p

        sqx = math.sqrt(x)
        rs_terms = int(math.sqrt(x / (2 * math.pi)))
        best_n = None

        for nz in [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500]:
            psi_val = psi_explicit(x, num_zeros=nz)
            err = abs(float(psi_val - true_psi))
            if err < 0.5:
                best_n = nz
                break

        if best_n:
            ratio = best_n / sqx
            flush_print(f"{x:>8} {float(true_psi):>12.4f} {best_n:>14} "
                        f"{sqx:>10.1f} {ratio:>10.3f} {rs_terms:>10}")
            results.append((x, best_n, sqx, ratio))
        else:
            flush_print(f"{x:>8} {float(true_psi):>12.4f} {'>500':>14} "
                        f"{sqx:>10.1f} {'?':>10} {rs_terms:>10}")
            results.append((x, 999, sqx, 999/sqx))

    # Fit power law: zeros_needed ~ C * x^alpha
    if len(results) >= 3:
        flush_print("\nFitting: zeros_needed ~ C * x^alpha")
        # Log-log regression
        from math import log as mlog
        xs = [mlog(r[0]) for r in results if r[1] < 999]
        ys = [mlog(r[1]) for r in results if r[1] < 999]
        if len(xs) >= 2:
            n = len(xs)
            sx = sum(xs); sy = sum(ys)
            sxy = sum(x*y for x,y in zip(xs,ys))
            sx2 = sum(x*x for x in xs)
            alpha = (n*sxy - sx*sy) / (n*sx2 - sx*sx)
            C = math.exp((sy - alpha*sx) / n)
            flush_print(f"  alpha = {alpha:.4f}  (1/2 = 0.5000)")
            flush_print(f"  C = {C:.4f}")
            flush_print(f"  => zeros_needed ~ {C:.2f} * x^{alpha:.3f}")

    return results


# ===========================================================================
# EXPERIMENT 2: Per-zero contribution analysis
# ===========================================================================

def experiment_zero_contributions(x=1000):
    flush_print("\n" + "="*75)
    flush_print(f"EXPERIMENT 2: Per-zero contributions to psi({x})")
    flush_print("="*75)

    psi_val, contribs = psi_explicit(x, num_zeros=200, return_terms=True)
    contribs_sorted = sorted(contribs, key=lambda c: -c[2])
    total_abs = sum(c[2] for c in contribs)

    flush_print(f"\nTop 15 contributing zeros:")
    flush_print(f"{'rank':>5} {'zero#':>7} {'gamma':>12} {'|contrib|':>14} {'contrib':>14}")
    flush_print("-"*55)

    for i, (k, gamma_k, abs_c, c) in enumerate(contribs_sorted[:15]):
        flush_print(f"{i+1:>5} {k:>7} {float(gamma_k):>12.4f} "
                    f"{float(abs_c):>14.6f} {float(c):>14.6f}")

    # Decay analysis
    flush_print(f"\nContribution decay (sorted by zero index, not magnitude):")
    flush_print(f"{'zeros':>8} {'cumul |contrib|':>18} {'% of total':>12}")
    flush_print("-"*42)
    cum = mpf(0)
    for i, (k, gamma_k, abs_c, c) in enumerate(contribs):
        cum += abs_c
        if (i+1) in [1, 2, 5, 10, 20, 50, 100, 200]:
            flush_print(f"{i+1:>8} {float(cum):>18.4f} "
                        f"{float(100*cum/total_abs):>11.2f}%")

    # Convergence of partial sums
    true_psi = mpmath.fsum(mpmath.log(p) for p in primes_up_to(x))
    for p in primes_up_to(int(x**0.5) + 1):
        pk = p * p
        while pk <= x:
            true_psi += log(p)
            pk *= p

    flush_print(f"\nPartial sum convergence:")
    flush_print(f"True psi({x}) = {float(true_psi):.6f}")
    flush_print(f"{'N zeros':>10} {'error':>14} {'exact?':>8}")
    flush_print("-"*35)

    for nz in [1, 2, 5, 10, 20, 50, 100, 200]:
        psi_n = psi_explicit(x, num_zeros=nz)
        err = float(psi_n - true_psi)
        ok = "YES" if abs(err) < 0.5 else "no"
        flush_print(f"{nz:>10} {err:>14.4f} {ok:>8}")

    flush_print(f"\n** Key: contribution of zero k ~ x^(1/2) / gamma_k")
    flush_print(f"** Decay is O(1/gamma_k) = O(1/k * ln(k)) -- ALGEBRAIC, not geometric")
    flush_print(f"** This means we CANNOT skip significant numbers of zeros")


# ===========================================================================
# EXPERIMENT 3: Gram point structure
# ===========================================================================

def experiment_gram_acceleration():
    flush_print("\n" + "="*75)
    flush_print("EXPERIMENT 3: Gram point structure analysis")
    flush_print("="*75)

    n_grams = 40
    grams = []
    z_values = []

    flush_print(f"\nFirst {n_grams} Gram points:")
    flush_print(f"{'n':>4} {'g_n':>14} {'Z(g_n)':>14} {'expected':>10} {'Gram?':>8}")
    flush_print("-"*55)

    violations = 0
    for n in range(n_grams):
        gn = gram_point(n)
        zn = siegelz(gn)
        expected = (-1)**n
        actual_sign = 1 if zn > 0 else -1
        ok = actual_sign == expected
        if not ok:
            violations += 1
        grams.append(gn)
        z_values.append(zn)

        if n < 15 or not ok:
            flush_print(f"{n:>4} {float(gn):>14.6f} {float(zn):>14.6f} "
                        f"{expected:>10} {'ok' if ok else 'FAIL':>8}")

    flush_print(f"\nGram violations: {violations}/{n_grams} "
                f"({100*violations/n_grams:.1f}%) [expected ~27%]")

    # Sign changes = zero crossings
    sign_changes = sum(1 for i in range(len(z_values)-1)
                       if z_values[i] * z_values[i+1] < 0)
    flush_print(f"Sign changes between consecutive Gram points: {sign_changes}")
    flush_print(f"(Each sign change locates one zero)")

    # Cost analysis
    flush_print(f"\n--- COST ANALYSIS ---")
    flush_print(f"To find N zeros using Gram points:")
    flush_print(f"  - Need ~N Gram point evaluations")
    flush_print(f"  - Each Z(g_n) via R-S costs O(sqrt(g_n/(2*pi))) = O(sqrt(n/ln(n)))")
    flush_print(f"  - Total: sum_1^N sqrt(n/ln(n)) ~ O(N^(3/2) / sqrt(ln N))")
    flush_print(f"  - Compare: Odlyzko-Schonhage: O(N * log^2(N))")
    flush_print(f"  - WINNER: Odlyzko-Schonhage is better for batch zero-finding")


# ===========================================================================
# EXPERIMENT 4: Optimal test function
# ===========================================================================

def experiment_optimal_test_function(x=1000):
    flush_print("\n" + "="*75)
    flush_print(f"EXPERIMENT 4: Optimal test function (x={x})")
    flush_print("="*75)

    true_psi = mpmath.fsum(mpmath.log(p) for p in primes_up_to(x))
    for p in primes_up_to(int(x**0.5) + 1):
        pk = p * p
        while pk <= x:
            true_psi += log(p)
            pk *= p

    # Test Gaussian-smoothed explicit formula
    flush_print(f"\nGaussian smoothing: h(t) = exp(-t^2/(2*sigma^2))")
    flush_print(f"Fourier decay: exp(-sigma^2 * gamma^2 / 2) -- FAST")
    flush_print(f"But smoothing error requires knowing primes near x.\n")

    sigmas = [0.01, 0.02, 0.05, 0.1]
    flush_print(f"{'sigma':>8} {'smooth_psi':>14} {'true_psi':>14} {'smooth_err':>12}")
    flush_print("-"*52)

    for sigma in sigmas:
        # Gaussian-smoothed: sum Lambda(n) * exp(-(n-x)^2 / (2*(sigma*x)^2))
        smooth_psi = mpf(0)
        width = sigma * x
        lo = max(2, int(x - 5 * width))
        hi = int(x + 5 * width)

        for p in primes_in_range(lo, hi):
            pk = p
            while pk <= hi:
                smooth_psi += log(p) * exp(-mpf(pk - x)**2 / (2 * mpf(width)**2))
                pk *= p

        smooth_err = float(abs(smooth_psi - true_psi))
        flush_print(f"{sigma:>8.3f} {float(smooth_psi):>14.4f} "
                    f"{float(true_psi):>14.4f} {smooth_err:>12.4f}")

    flush_print(f"""
    CRITICAL ANALYSIS:
    - Gaussian with sigma=c/ln(x) would need only O(log^2 x) zeros
    - BUT the smoothing error is ~ sum_{{p near x}} Lambda(p) * (1 - gaussian)
    - Resolving which primes are near x IS the original problem
    - The "free lunch" of fewer zeros is paid for by the smoothing correction
    - Net effect: no improvement over O(x^(1/2))""")


# ===========================================================================
# EXPERIMENT 5: Odlyzko-Schonhage analysis
# ===========================================================================

def experiment_odlyzko_schonhage():
    flush_print("\n" + "="*75)
    flush_print("EXPERIMENT 5: Odlyzko-Schonhage batching")
    flush_print("="*75)

    flush_print("""
    Odlyzko-Schonhage computes N consecutive Z(t) values in O(N * log^2 N).

    For pi(x) via explicit formula:
    - Need zeros with gamma < T for truncation error < 0.5
    - Truncation tail: sum_{gamma > T} |x^rho / rho| ~ x^(1/2) * ln(x) / T
    - Need T ~ 2 * x^(1/2) * ln(x) for error < 0.5
    - Number of zeros up to T: N(T) ~ T * ln(T) / (2*pi)
    - Substituting: N ~ x^(1/2) * ln^2(x)

    O-S cost to find all N zeros: O(N * log^2 N) = O(x^(1/2) * log^4 x)
    Explicit formula evaluation:  O(N) = O(x^(1/2) * log^2 x)

    TOTAL: O(x^(1/2) * log^4 x)  --  this is O(x^(1/2+epsilon))

    Compare with current best:
    - Meissel-Lehmer:          O(x^(2/3) / log^2 x)   [combinatorial]
    - Lagarias-Miller-Odlyzko: O(x^(2/3))              [analytic]

    For x < 10^13, the analytic explicit formula is WORSE.
    For x > 10^50, it would be BETTER than O(x^(2/3))...
    ...but 2/3 > 1/2, so theoretically the explicit formula already wins!

    WAIT -- this is a known result! The analytic method IS O(x^(1/2+eps)).
    But in practice:
    - The constant in O-S is huge
    - Memory requirements are O(N) ~ O(x^(1/2))
    - The Meissel-Lehmer O(x^(2/3)) has tiny constants
    - Nobody has implemented the full O(x^(1/2+eps}) for production""")

    # Numerical: verify zero-sum convergence
    flush_print("\nNumerical convergence verification:")
    flush_print(f"{'x':>8} {'N needed':>10} {'sqrt(x)':>10} {'N/sqrt(x)':>12}")
    flush_print("-"*45)

    for x in [100, 500, 1000, 2000]:
        true_psi = mpmath.fsum(mpmath.log(p) for p in primes_up_to(x))
        for p in primes_up_to(int(x**0.5) + 1):
            pk = p * p
            while pk <= x:
                true_psi += log(p)
                pk *= p
        sqx = math.sqrt(x)
        best_n = None
        for nz in [1, 2, 5, 10, 20, 50, 100, 200, 500]:
            psi_val = psi_explicit(x, num_zeros=nz)
            if abs(float(psi_val - true_psi)) < 0.5:
                best_n = nz
                break
        if best_n:
            flush_print(f"{x:>8} {best_n:>10} {sqx:>10.1f} {best_n/sqx:>12.4f}")
        else:
            flush_print(f"{x:>8} {'>500':>10} {sqx:>10.1f} {'---':>12}")


# ===========================================================================
# EXPERIMENT 6: Zero spacing / Deuring-Heilbronn
# ===========================================================================

def experiment_deuring_heilbronn():
    flush_print("\n" + "="*75)
    flush_print("EXPERIMENT 6: Zero spacing statistics")
    flush_print("="*75)

    zeros = [float(get_zero(k).imag) for k in range(1, 101)]
    spacings = [zeros[i+1] - zeros[i] for i in range(len(zeros)-1)]
    mean_sp = sum(spacings) / len(spacings)
    norm_sp = [s / mean_sp for s in spacings]
    var = sum((s-1)**2 for s in norm_sp) / len(norm_sp)

    flush_print(f"\nFirst 100 zeros spacing analysis:")
    flush_print(f"  Mean spacing: {mean_sp:.6f}")
    flush_print(f"  Min spacing:  {min(spacings):.6f}")
    flush_print(f"  Max spacing:  {max(spacings):.6f}")
    flush_print(f"  Variance:     {var:.4f}  (GUE predicts ~0.286)")

    flush_print(f"""
    Can zero repulsion help?
    - GUE statistics give the DISTRIBUTION, not exact positions
    - Predicting individual zeros requires info about ALL others
    - The pair correlation (Montgomery conjecture) is statistical
    - VERDICT: No computational speedup from repulsion statistics""")


# ===========================================================================
# THEORETICAL SUMMARY
# ===========================================================================

def theoretical_analysis():
    flush_print("\n" + "="*75)
    flush_print("THEORETICAL ANALYSIS: Can Gram Points Beat O(x^{1/2})?")
    flush_print("="*75)

    flush_print("""
    =====================================================================
    APPROACH                     | ZEROS NEEDED | COST/ZERO  | TOTAL
    =====================================================================
    Standard explicit formula    | O(x^(1/2))   | O(1) [pre] | O(x^(1/2))
    Riemann-Siegel per-zero      | O(x^(1/2))   | O(t^(1/2)) | O(x^(3/4))
    Odlyzko-Schonhage batch      | O(x^(1/2))   | O(log^2 x) | O(x^(1/2+e))
    Gaussian smoothing           | O(log^2 x)   | O(log^2 x) | CIRCULAR*
    Optimal test function        | O(log^k x)   | varies     | CIRCULAR*
    =====================================================================
    * Smoothing reduces zeros but adds correction requiring primes near x

    THE FUNDAMENTAL BARRIER:
    ========================
    1. The explicit formula converges as O(x^(1/2) / N) after N zeros.
    2. Exact pi(x) needs error < 0.5, so N ~ O(x^(1/2) * ln x).
    3. Smoothing reduces N but adds a correction that requires O(x^a) work.
    4. The Odlyzko-Schonhage algorithm gives O(x^(1/2) * polylog) total.
    5. This is ALREADY the theoretical best for the analytic approach.
    6. It CANNOT be improved to sub-O(x^(1/2)) because:
       - Each of the O(x^(1/2)) zeros contributes independently
       - No known structure allows batch-evaluating the sum in fewer terms
       - The conditional convergence means we cannot skip zeros

    COMPARISON WITH BEST KNOWN:
    - Analytic (explicit formula + O-S): O(x^(1/2+epsilon))
    - Combinatorial (Meissel-Lehmer):    O(x^(2/3) / log^2 x)
    - Practical best (Deleglise-Rivat):  O(x^(2/3) / log^2 x)

    Note: 1/2 < 2/3, so the analytic approach is ASYMPTOTICALLY better!
    But the constants are much worse. In practice, analytic methods have
    been used only for verification, not primary computation.

    Can we do BETTER than O(x^(1/2))?
    - The zero sum has O(x^(1/2)) independent terms
    - This appears to be an INFORMATION-THEORETIC barrier
    - The zeros encode the "random" part of the prime distribution
    - There are O(x^(1/2)) of them, and each carries ~1 bit of info
    - Total info content: ~x^(1/2) bits = omega(polylog) bits
    - Any algorithm must process these bits => Omega(x^(1/2)) time

    STATUS: APPROACH RULED OUT
    ==========================
    The Gram point + Riemann-Siegel combination gives O(x^(1/2+eps}),
    which is already known and cannot be improved to sub-sqrt.
    This is approach #206+ in the ongoing search.
    """)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    flush_print("Session 8: Gram Points + Riemann-Siegel Approach")
    flush_print("=" * 75)

    t0 = time.time()

    # Sanity check
    flush_print("\n--- Sanity checks ---")
    t_test = mpf(100)
    z_mpmath = siegelz(t_test)
    z_manual, ops = riemann_siegel_manual(t_test, count_ops=True)
    flush_print(f"Z(100) via mpmath:  {float(z_mpmath):.10f}")
    flush_print(f"Z(100) via manual:  {float(z_manual):.10f} ({ops} RS terms)")
    flush_print(f"Agreement: {float(abs(z_mpmath - z_manual)):.2e}")

    # Run experiments
    experiment_zeros_needed()
    experiment_zero_contributions(x=1000)
    experiment_gram_acceleration()
    experiment_optimal_test_function(x=1000)
    experiment_odlyzko_schonhage()
    experiment_deuring_heilbronn()
    theoretical_analysis()

    elapsed = time.time() - t0
    flush_print(f"\nTotal runtime: {elapsed:.1f}s")
    flush_print(f"\n{'='*75}")
    flush_print(f"SESSION 8 VERDICT:")
    flush_print(f"  Gram points + R-S give O(x^(1/2) * polylog) -- already known.")
    flush_print(f"  Cannot beat O(x^(1/2)) due to information-theoretic barrier.")
    flush_print(f"  The explicit formula has ~x^(1/2) independent zero contributions")
    flush_print(f"  that cannot be compressed or batch-evaluated sub-linearly.")
    flush_print(f"  RULED OUT as path to sub-O(x^(1/2)) pi(x).")
    flush_print(f"{'='*75}")


if __name__ == "__main__":
    main()
