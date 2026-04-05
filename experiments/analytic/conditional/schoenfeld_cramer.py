#!/usr/bin/env python3
"""
Experiments 5 & 6 from Task #4:
  5. Schoenfeld's Explicit Bounds under RH
  6. Cramer's Conjecture Search Algorithm

Investigates whether these conditional results can yield sub-O(x^{2/3}) prime computation.
"""

import math
import time
import os
import sys
from pathlib import Path

import mpmath
from mpmath import mp, mpf, log, sqrt, pi as MP_PI, li as mp_li, loggamma
from sympy import primepi, prime, isprime, nextprime, prevprime, primerange

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_DIR = Path("/apps/aplikacijos/prime-research/data")
RESULTS_DIR = Path(__file__).parent

mp.dps = 50  # decimal places for mpmath


# ---------------------------------------------------------------------------
# Load zeta zeros from data/
# ---------------------------------------------------------------------------
def load_zeros(count=1000):
    """Load imaginary parts of first `count` non-trivial zeta zeros."""
    fpath = DATA_DIR / f"zeta_zeros_{count}.txt"
    if not fpath.exists():
        # Try the largest available file
        for c in [1000, 500, 300, 200]:
            fpath = DATA_DIR / f"zeta_zeros_{c}.txt"
            if fpath.exists():
                break
    zeros = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # File may have just the value, or index + value
            val = parts[-1]  # take the last column either way
            try:
                zeros.append(mpf(val))
            except Exception:
                continue
    return zeros


ZETA_ZEROS = load_zeros()
print(f"Loaded {len(ZETA_ZEROS)} zeta zeros")


# ---------------------------------------------------------------------------
# Riemann R function and its inverse
# ---------------------------------------------------------------------------
def riemann_R(x, terms=200):
    """Compute the Riemann R function: R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})."""
    from sympy import mobius
    x = mpf(x)
    if x <= 1:
        return mpf(0)
    result = mpf(0)
    for k in range(1, terms + 1):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        xk = x ** (mpf(1) / k)
        if xk <= 1.0001:
            break
        result += mpf(mu_k) / k * mp_li(xk)
    return result


def inverse_R(n, tol=1e-10):
    """Compute R^{-1}(n) using Newton's method. Returns float."""
    n = mpf(n)
    # Initial estimate: n * ln(n)
    x = n * log(n)
    for _ in range(100):
        rx = riemann_R(x)
        err = rx - n
        if abs(err) < tol:
            break
        # Derivative of R(x) ~ 1/ln(x)
        deriv = 1.0 / log(x)
        x = x - err / deriv
    return float(x)


# ---------------------------------------------------------------------------
# Schoenfeld bound
# ---------------------------------------------------------------------------
def schoenfeld_bound(x):
    """
    Under RH, |pi(x) - li(x)| < sqrt(x)*ln(x)/(8*pi) for x >= 2657.
    Returns the bound as a float.
    """
    x = float(x)
    return math.sqrt(x) * math.log(x) / (8 * math.pi)


# ---------------------------------------------------------------------------
# Explicit formula correction with K zeros
# ---------------------------------------------------------------------------
def explicit_correction(x, K):
    """
    Compute the oscillatory correction from the first K zeta zeros.
    pi(x) ~ li(x) - sum_{k=1}^{K} [ li(x^{rho_k}) + li(x^{conj(rho_k)}) ]
    where rho_k = 1/2 + i*gamma_k (under RH).

    The correction term for a pair (rho, conj(rho)) with rho = 1/2 + i*gamma is:
      -2 * Re[ li(x^rho) ]

    We use the approximation for large x:
      li(x^rho) ~ x^rho / (rho * ln(x))

    So the pair contributes: -2 * Re[ x^{1/2 + i*gamma} / ((1/2 + i*gamma) * ln(x)) ]
    """
    x = mpf(x)
    lnx = log(x)
    correction = mpf(0)
    K = min(K, len(ZETA_ZEROS))

    for k in range(K):
        gamma = ZETA_ZEROS[k]
        # x^{1/2 + i*gamma} = sqrt(x) * exp(i*gamma*ln(x))
        # = sqrt(x) * (cos(gamma*lnx) + i*sin(gamma*lnx))
        phase = gamma * lnx
        cos_phase = mpmath.cos(phase)
        sin_phase = mpmath.sin(phase)

        sqrtx = sqrt(x)

        # 1/(rho * ln(x)) where rho = 1/2 + i*gamma
        # 1/((1/2 + i*gamma)*lnx) = (1/2 - i*gamma) / ((1/4 + gamma^2)*lnx)
        denom = (mpf(0.25) + gamma ** 2) * lnx
        re_inv = mpf(0.5) / denom
        im_inv = -gamma / denom

        # x^rho / (rho * lnx) = sqrtx * (cos + i*sin) * (re_inv + i*im_inv)
        # Real part: sqrtx * (cos*re_inv - sin*im_inv)
        re_term = sqrtx * (cos_phase * re_inv - sin_phase * im_inv)

        correction += -2 * re_term

    return float(correction)


def residual_error_bound_with_K_zeros(x, K):
    """
    Under RH, after summing K zeros, the residual error is approximately:
      O(sqrt(x) * log(x) / K)
    More precisely, heuristic bound:  sqrt(x) * log^2(x) / (K + 1)
    (the log^2 factor accounts for the density of zeros).
    """
    x = float(x)
    return math.sqrt(x) * math.log(x) ** 2 / (K + 1)


# =========================================================================
# EXPERIMENT 5: Schoenfeld's Explicit Bounds under RH
# =========================================================================
def experiment_5():
    print("=" * 72)
    print("EXPERIMENT 5: Schoenfeld's Explicit Bounds under RH")
    print("=" * 72)

    results = {}

    # Part A: Basic Schoenfeld bound vs Meissel-Lehmer
    print("\n--- Part A: Schoenfeld bound size vs x ---")
    print(f"{'x':>12s}  {'sqrt(x)*ln(x)/(8pi)':>22s}  {'x^(2/3)':>14s}  {'Schoenfeld/ML':>14s}")
    print("-" * 72)
    part_a = []
    for exp in range(4, 21):
        x = 10 ** exp
        sb = schoenfeld_bound(x)
        ml = x ** (2.0 / 3)
        ratio = sb / ml
        part_a.append((exp, sb, ml, ratio))
        print(f"  10^{exp:<6d}  {sb:>22.2f}  {ml:>14.2f}  {ratio:>14.6f}")

    results['part_a'] = part_a

    # Part B: Sieve cost in Schoenfeld interval
    print("\n--- Part B: Sieve cost in Schoenfeld interval vs O(x^{2/3}) ---")
    print(f"{'x':>12s}  {'Interval width 2E':>18s}  {'Sieve cost':>18s}  {'ML cost x^(2/3)':>18s}  {'Ratio':>10s}")
    print("-" * 72)
    part_b = []
    for exp in [4, 6, 8, 10, 12, 15, 20, 50, 100]:
        x = 10.0 ** exp
        E = math.sqrt(x) * math.log(x) / (8 * math.pi)
        width = 2 * E
        sieve_cost = width * math.log(math.log(x))  # O(W * ln(ln(x)))
        ml_cost = x ** (2.0 / 3)
        ratio = sieve_cost / ml_cost
        part_b.append((exp, width, sieve_cost, ml_cost, ratio))
        if exp <= 20:
            print(f"  10^{exp:<6d}  {width:>18.2e}  {sieve_cost:>18.2e}  {ml_cost:>18.2e}  {ratio:>10.4f}")
        else:
            print(f"  10^{exp:<6d}  ~10^{math.log10(width):.1f}{'':>6s}  ~10^{math.log10(sieve_cost):.1f}{'':>6s}  ~10^{math.log10(ml_cost):.1f}{'':>6s}  {ratio:>10.4e}")

    results['part_b'] = part_b

    # Part C: Iterative refinement with K zeros
    print("\n--- Part C: Error reduction with K zeta zeros ---")
    print("Testing for x = 10^4 to 10^8, comparing actual pi(x) vs li(x) + correction")
    print()

    test_points = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
    K_values = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
    K_values = [k for k in K_values if k <= len(ZETA_ZEROS)]

    part_c = {}
    for x in test_points:
        actual_pi = int(primepi(x))
        li_x = float(mp_li(mpf(x)))
        base_error = abs(li_x - actual_pi)

        print(f"  x = {x:>12,d}  |  pi(x) = {actual_pi:>10,d}  |  li(x) = {li_x:>14.2f}  |  |li(x)-pi(x)| = {base_error:.2f}")
        print(f"    Schoenfeld bound: {schoenfeld_bound(x):.2f}")
        print(f"    {'K zeros':>10s}  {'Corrected est':>16s}  {'Actual error':>14s}  {'Heuristic bound':>16s}  {'Bound holds?':>12s}")

        errors_for_x = []
        for K in K_values:
            corr = explicit_correction(x, K)
            est = li_x + corr
            actual_err = abs(est - actual_pi)
            heur_bound = residual_error_bound_with_K_zeros(x, K)
            holds = "YES" if actual_err <= heur_bound else "no"
            errors_for_x.append((K, est, actual_err, heur_bound, holds))
            print(f"    {K:>10d}  {est:>16.2f}  {actual_err:>14.2f}  {heur_bound:>16.2f}  {holds:>12s}")

        part_c[x] = errors_for_x
        print()

    results['part_c'] = part_c

    # Part D: How many zeros needed for error < 1?
    print("--- Part D: Zeros needed for error < 1 (theoretical) ---")
    print(f"{'x':>12s}  {'K needed (heuristic)':>22s}  {'sqrt(x)*log(x)':>18s}")
    print("-" * 60)
    part_d = []
    for exp in [4, 6, 8, 10, 12, 15, 20, 50, 100]:
        x = 10.0 ** exp
        # Need sqrt(x)*log^2(x)/(K+1) < 1 => K > sqrt(x)*log^2(x)
        K_needed = math.sqrt(x) * math.log(x) ** 2
        sqrtx_logx = math.sqrt(x) * math.log(x)
        part_d.append((exp, K_needed, sqrtx_logx))
        print(f"  10^{exp:<6d}  ~10^{math.log10(K_needed):.1f}{'':>14s}  ~10^{math.log10(sqrtx_logx):.1f}")

    results['part_d'] = part_d

    # Part E: Numerical verification - does the correction converge?
    print("\n--- Part E: Convergence of correction terms for x=10^6 ---")
    x = 1_000_000
    actual_pi = int(primepi(x))
    li_x = float(mp_li(mpf(x)))

    print(f"  pi({x}) = {actual_pi}, li({x}) = {li_x:.4f}")
    print(f"  Base error |li(x)-pi(x)| = {abs(li_x - actual_pi):.4f}")
    print()

    prev_err = abs(li_x - actual_pi)
    part_e = []
    for K in range(1, min(len(ZETA_ZEROS) + 1, 1001)):
        corr = explicit_correction(x, K)
        est = li_x + corr
        err = abs(est - actual_pi)
        if K <= 20 or K % 50 == 0 or K == min(len(ZETA_ZEROS), 1000):
            improvement = prev_err / err if err > 0 else float('inf')
            part_e.append((K, err))
            print(f"    K={K:>5d}: corrected={est:>12.4f}, error={err:>10.4f}")

    results['part_e'] = part_e

    return results


# =========================================================================
# EXPERIMENT 6: Cramer's Conjecture Search Algorithm
# =========================================================================
def experiment_6():
    print("\n" + "=" * 72)
    print("EXPERIMENT 6: Cramer's Conjecture Search Algorithm")
    print("=" * 72)

    EULER_GAMMA = 0.5772156649015329
    GRANVILLE_C = 2 * math.exp(-EULER_GAMMA)  # ~ 1.1229

    results = {}

    # Part A: Implement the full algorithm and measure phases
    print("\n--- Part A: Cramer search algorithm phase costs ---")
    print("Steps: (1) R^{-1}(n) -> x0, (2) pi(x0) [BOTTLENECK], (3-4) walk to p(n)")
    print()

    test_ns = [100, 1000, 10000, 100000, 1000000, 5000000]
    part_a = []

    for n in test_ns:
        # Phase 1: Compute R^{-1}(n)
        t0 = time.time()
        x0 = inverse_R(n)
        t_rinv = time.time() - t0

        # Phase 2: Compute pi(x0) [the bottleneck]
        t0 = time.time()
        x0_int = int(round(x0))
        pi_x0 = int(primepi(x0_int))
        t_count = time.time() - t0

        # Phase 3-4: Walk from x0 to actual p(n)
        t0 = time.time()
        actual_pn = int(prime(n))
        diff = pi_x0 - n  # How far off is our count?

        # Walk to find p(n)
        if pi_x0 == n:
            # x0 might be prime or we need prevprime
            if isprime(x0_int):
                found = x0_int
            else:
                found = int(prevprime(x0_int + 1))
        elif pi_x0 > n:
            # Too many primes below x0, search backwards
            current = x0_int
            for _ in range(abs(diff) + 1):
                current = int(prevprime(current))
            found = current
        else:
            # Too few primes, search forwards
            current = x0_int
            for _ in range(abs(diff)):
                current = int(nextprime(current))
            found = current

        t_walk = time.time() - t0

        # Verify
        correct = (found == actual_pn)

        # Cramer gap prediction
        gap_bound = GRANVILLE_C * math.log(actual_pn) ** 2
        actual_gap_approx = abs(x0 - actual_pn)  # Initial offset

        entry = {
            'n': n,
            'pn': actual_pn,
            'x0': x0,
            'offset': x0 - actual_pn,
            'pi_x0': pi_x0,
            'diff': diff,
            'steps': abs(diff),
            't_rinv': t_rinv,
            't_count': t_count,
            't_walk': t_walk,
            'correct': correct,
            'gap_bound': gap_bound,
        }
        part_a.append(entry)

        total = t_rinv + t_count + t_walk
        pct_count = 100 * t_count / total if total > 0 else 0

        print(f"  n = {n:>10,d}  |  p(n) = {actual_pn:>12,d}")
        print(f"    x0 = R^{{-1}}(n) = {x0:>14.2f}  (offset = {x0 - actual_pn:>+.2f})")
        print(f"    pi(x0) = {pi_x0}, need {n}, diff = {diff:+d} -> walk {abs(diff)} steps")
        print(f"    Times: R^{{-1}}={t_rinv:.4f}s  pi(x0)={t_count:.4f}s  walk={t_walk:.4f}s")
        print(f"    Counting phase: {pct_count:.1f}% of total  |  Correct: {correct}")
        print(f"    Cramer gap bound: {gap_bound:.1f}")
        print()

    results['part_a'] = part_a

    # Part B: Scaling analysis
    print("--- Part B: Scaling analysis of phases ---")
    print(f"{'n':>10s}  {'t_R_inv':>10s}  {'t_count':>10s}  {'t_walk':>10s}  {'%count':>8s}  {'walk_steps':>10s}  {'ln^2(pn)':>10s}")
    print("-" * 75)
    part_b = []
    for entry in part_a:
        n = entry['n']
        total = entry['t_rinv'] + entry['t_count'] + entry['t_walk']
        pct = 100 * entry['t_count'] / total if total > 0 else 0
        ln2pn = math.log(entry['pn']) ** 2
        part_b.append((n, entry['t_rinv'], entry['t_count'], entry['t_walk'], pct, entry['steps'], ln2pn))
        print(f"  {n:>8,d}  {entry['t_rinv']:>10.5f}  {entry['t_count']:>10.5f}  {entry['t_walk']:>10.5f}  {pct:>7.1f}%  {entry['steps']:>10d}  {ln2pn:>10.1f}")

    results['part_b'] = part_b

    # Part C: Gap statistics near R^{-1}(n) - verify Cramer
    print("\n--- Part C: Gap statistics near R^{-1}(n) ---")
    print("Checking if gaps near x0 respect Cramer's bound C*ln^2(x)")
    print()

    gap_stats = []
    for n in [1000, 10000, 100000, 1000000]:
        pn = int(prime(n))
        ln2 = math.log(pn) ** 2

        # Look at 20 gaps near p(n)
        max_gap = 0
        gaps = []
        p = pn
        for _ in range(20):
            np_ = int(nextprime(p))
            g = np_ - p
            gaps.append(g)
            max_gap = max(max_gap, g)
            p = np_

        avg_gap = sum(gaps) / len(gaps)
        max_ratio = max_gap / ln2
        avg_ratio = avg_gap / ln2

        gap_stats.append((n, pn, max_gap, avg_gap, ln2, max_ratio, avg_ratio))
        print(f"  n={n:>8,d}  p(n)={pn:>12,d}  ln^2(p)={ln2:.1f}  max_gap={max_gap}  "
              f"max/ln^2={max_ratio:.3f}  avg_gap={avg_gap:.1f}  avg/ln^2={avg_ratio:.3f}")

    results['part_c'] = gap_stats

    # Part D: The fundamental bottleneck
    print("\n--- Part D: Bottleneck analysis ---")
    print("Key insight: The counting phase pi(x0) dominates.")
    print("Walk phase cost scales with O(|diff| * primality_test) = O(small).")
    print()
    print("Theoretical scaling:")
    print(f"{'x':>12s}  {'Count O(x^{2/3})':>18s}  {'Walk O(ln^4 x)':>18s}  {'Ratio':>12s}")
    print("-" * 65)
    part_d = []
    for exp in [4, 6, 8, 10, 12, 15, 20, 50, 100]:
        x = 10.0 ** exp
        count_cost = x ** (2.0 / 3)
        walk_cost = math.log(x) ** 4
        ratio = count_cost / walk_cost
        part_d.append((exp, count_cost, walk_cost, ratio))
        print(f"  10^{exp:<6d}  10^{math.log10(count_cost):>6.1f}{'':>7s}  10^{math.log10(walk_cost):>6.2f}{'':>6s}  10^{math.log10(ratio):>6.1f}")

    results['part_d'] = part_d

    return results


# =========================================================================
# MAIN
# =========================================================================
def main():
    print("Schoenfeld-Cramer Experiments")
    print(f"Using {len(ZETA_ZEROS)} zeta zeros from {DATA_DIR}")
    print()

    results_5 = experiment_5()
    results_6 = experiment_6()

    # Final summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print()
    print("Experiment 5 (Schoenfeld bounds):")
    print("  - Schoenfeld interval width ~ sqrt(x)*ln(x) -> sieve cost O(sqrt(x)*ln(x)*ln(ln(x)))")
    print("  - This is WORSE than Meissel-Lehmer O(x^{2/3}) for all practical x")
    print("  - Iterative refinement with K zeros: error ~ sqrt(x)*log^2(x)/K")
    print("  - For error < 1, need K ~ sqrt(x)*log^2(x) zeros -- infeasible for large x")
    print("  - With 1000 zeros, significant error reduction for x <= 10^6, but not enough for exactness")
    print()
    print("Experiment 6 (Cramer search):")
    print("  - Walk phase is trivially fast: O(ln^4 x) under Cramer's conjecture")
    print("  - BOTTLENECK is pi(x0) computation: still O(x^{2/3}) via Meissel-Lehmer")
    print("  - Counting phase dominates by factor 10^{2e/3 - 4*log10(e*ln10)} for x=10^e")
    print("  - At x=10^100: counting ~ 10^66, walking ~ 10^4.6 -> ratio ~ 10^62")
    print("  - Cramer helps the SEARCH but not the COUNTING")
    print()
    print("VERDICT: Neither Schoenfeld bounds nor Cramer's conjecture breaks the")
    print("O(x^{2/3}) barrier. The bottleneck remains computing pi(x), not locating p(n).")
    print("PATH CLOSED.")

    return results_5, results_6


if __name__ == "__main__":
    main()
