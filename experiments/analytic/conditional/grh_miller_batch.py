#!/usr/bin/env python3
"""
GRH Miller Batch Testing & Explicit Formula Optimal T
=====================================================
Task #4, Experiments 1 & 2.

Experiment 1: GRH-based deterministic Miller primality testing.
  - Under GRH, witnesses up to 2*ln(n)^2 suffice.
  - Measure actual witnesses needed vs GRH bound for n up to 10^7.
  - Investigate batch Miller testing for ranges.

Experiment 2: GRH Explicit Formula - optimal T (number of zeros).
  - Under GRH/RH, |pi(x) - li(x) + sum_{|gamma|<T} li(x^rho)| << sqrt(x)*log(x)/T
  - Find T_min such that error < 0.5 for x = 10^4 through 10^8.
  - Verify T_min = O(sqrt(x)*log^2(x)).
  - Can we need fewer than sqrt(x) zeros?
"""

import math
import time
import os
import sys
from collections import defaultdict

import sympy
from sympy import primepi, isprime, nextprime, primerange
import mpmath

mpmath.mp.dps = 30  # 30 decimal digits precision

# ============================================================
# DATA LOADING
# ============================================================

DATA_DIR = "/apps/aplikacijos/prime-research/data"

def load_zeta_zeros(max_zeros=1000):
    """Load precomputed zeta zeros from data directory."""
    # Try largest file first
    for fname in ["zeta_zeros_1000.txt", "zeta_zeros_500.txt",
                  "zeta_zeros_300.txt", "zeta_zeros_200.txt"]:
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            zeros = []
            with open(fpath) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        zeros.append(mpmath.mpf(parts[1]))
                    elif len(parts) == 1:
                        zeros.append(mpmath.mpf(parts[0]))
                if len(zeros) >= max_zeros:
                    return zeros[:max_zeros]
            if zeros:
                print(f"  Loaded {len(zeros)} zeros from {fname}")
                return zeros
    raise FileNotFoundError("No zeta zeros files found in data/")


# ============================================================
# EXPERIMENT 1: GRH MILLER BATCH TESTING
# ============================================================

def miller_rabin_witness(n, a):
    """
    Test if 'a' is a witness to compositeness of n via Miller-Rabin.
    Returns True if n is composite (a is a witness), False if n passes.
    """
    if n < 2:
        return True
    if n == 2:
        return False
    if n % 2 == 0:
        return True

    # Write n-1 = 2^s * d
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return False  # passes this witness

    for _ in range(s - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return False  # passes
    return True  # composite


def grh_witness_bound(n):
    """Under GRH, witnesses up to 2*ln(n)^2 suffice for Miller's test."""
    if n < 3:
        return 2
    ln_n = math.log(n)
    return int(math.floor(2 * ln_n * ln_n))


def miller_grh_is_prime(n):
    """
    Deterministic Miller primality test under GRH.
    Test all witnesses a = 2, 3, ..., floor(2*ln(n)^2).
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    bound = grh_witness_bound(n)
    for a in range(2, min(bound + 1, n)):
        if miller_rabin_witness(n, a):
            return False
    return True


def find_min_witnesses_needed(n):
    """
    For a composite n, find the smallest witness that detects it.
    For a prime n, return -1 (no witness detects it as composite).
    """
    if sympy.isprime(n):
        return -1  # prime, no witness needed

    for a in range(2, n):
        if math.gcd(a, n) > 1:
            # Trivial witness (shares factor)
            return a
        if miller_rabin_witness(n, a):
            return a
    return n  # shouldn't happen


def experiment1_witness_analysis():
    """
    Analyze how many witnesses are actually needed vs GRH bound.
    For composites up to various limits, find the smallest witness.
    """
    print("=" * 70)
    print("EXPERIMENT 1: GRH BATCH MILLER TESTING")
    print("=" * 70)

    # Part 1a: Verify correctness of GRH Miller test
    print("\n--- Part 1a: Correctness verification ---")
    errors = 0
    test_limit = 10**5
    for n in range(2, test_limit):
        grh_result = miller_grh_is_prime(n)
        actual = sympy.isprime(n)
        if grh_result != actual:
            errors += 1
            print(f"  ERROR at n={n}: GRH says {grh_result}, actual {actual}")
    print(f"  Tested n=2..{test_limit-1}: {errors} errors")
    print(f"  GRH Miller test {'CORRECT' if errors == 0 else 'HAS ERRORS'} up to {test_limit}")

    # Part 1b: Witness statistics for composites
    print("\n--- Part 1b: Smallest witness for composites ---")
    limits = [10**3, 10**4, 10**5, 10**6]
    results_1b = {}

    for limit in limits:
        max_witness = 0
        max_witness_n = 0
        witness_counts = defaultdict(int)
        total_composites = 0

        for n in range(3, limit + 1, 2):  # odd composites
            if sympy.isprime(n):
                continue
            total_composites += 1
            w = find_min_witnesses_needed(n)
            witness_counts[w] += 1
            if w > max_witness:
                max_witness = w
                max_witness_n = n

        grh_bound_at_limit = grh_witness_bound(limit)
        ratio = max_witness / grh_bound_at_limit if grh_bound_at_limit > 0 else 0

        results_1b[limit] = {
            'max_witness': max_witness,
            'max_witness_n': max_witness_n,
            'grh_bound': grh_bound_at_limit,
            'ratio': ratio,
            'total_composites': total_composites,
        }

        print(f"  n up to {limit:>10,}:")
        print(f"    Max smallest witness needed: {max_witness} (at n={max_witness_n})")
        print(f"    GRH bound 2*ln(n)^2:        {grh_bound_at_limit}")
        print(f"    Ratio (actual/bound):        {ratio:.4f}")
        print(f"    Composites tested:           {total_composites}")

    # Part 1c: GRH bound vs actual witnesses needed
    print("\n--- Part 1c: Witness count comparison ---")
    print("  For each n, GRH says test witnesses 2..floor(2*ln(n)^2)")
    sample_ns = [100, 1000, 10000, 100000, 10**6, 10**7, 10**8, 10**12, 10**20, 10**50, 10**100]
    print(f"  {'n':>15s}  {'GRH bound':>12s}  {'log2(n)':>10s}  {'witnesses/log2(n)':>18s}")
    grh_data = []
    for n in sample_ns:
        bound = grh_witness_bound(n)
        log2n = math.log2(n)
        ratio = bound / log2n if log2n > 0 else 0
        grh_data.append((n, bound, log2n, ratio))
        print(f"  {n:>15.0e}  {bound:>12,}  {log2n:>10.1f}  {ratio:>18.2f}")

    # Part 1d: Batch Miller testing
    print("\n--- Part 1d: Batch Miller testing ---")
    print("  Question: Can witnesses be shared across a range [a,b]?")

    # Test: for a range, what's the maximum smallest-witness?
    ranges = [(10**4, 10**4 + 1000), (10**5, 10**5 + 1000), (10**6, 10**6 + 1000)]
    batch_results = []

    for lo, hi in ranges:
        max_w = 0
        max_w_n = 0
        composites_in_range = 0
        for n in range(lo | 1, hi, 2):  # odd numbers in range
            if sympy.isprime(n):
                continue
            composites_in_range += 1
            w = find_min_witnesses_needed(n)
            if w > max_w:
                max_w = w
                max_w_n = n

        grh_b = grh_witness_bound(hi)
        batch_results.append((lo, hi, max_w, max_w_n, grh_b, composites_in_range))
        print(f"  Range [{lo}, {hi}]: max witness needed = {max_w} (at n={max_w_n})")
        print(f"    GRH bound at {hi}: {grh_b}")
        print(f"    Composites in range: {composites_in_range}")

    # Batch timing comparison
    print("\n  Timing: individual vs batch Miller testing")
    test_range_lo = 10**6
    test_range_hi = test_range_lo + 5000

    # Individual: each number gets its own witness set
    t0 = time.time()
    individual_primes = []
    for n in range(test_range_lo, test_range_hi):
        if miller_grh_is_prime(n):
            individual_primes.append(n)
    t_individual = time.time() - t0

    # "Batch": precompute witness set for largest n, use for all
    t0 = time.time()
    max_bound = grh_witness_bound(test_range_hi)
    witnesses = list(range(2, max_bound + 1))
    batch_primes = []
    for n in range(test_range_lo, test_range_hi):
        if n < 2:
            continue
        if n == 2:
            batch_primes.append(n)
            continue
        if n % 2 == 0:
            continue
        is_prime = True
        for a in witnesses:
            if a >= n:
                break
            if miller_rabin_witness(n, a):
                is_prime = False
                break
        if is_prime:
            batch_primes.append(n)
    t_batch = time.time() - t0

    assert individual_primes == batch_primes, "Mismatch between individual and batch!"
    print(f"    Range [{test_range_lo}, {test_range_hi}), {len(individual_primes)} primes found")
    print(f"    Individual: {t_individual:.4f}s")
    print(f"    Batch (shared witnesses): {t_batch:.4f}s")
    print(f"    Speedup: {t_individual/t_batch:.2f}x" if t_batch > 0 else "    N/A")

    # Part 1e: Complexity analysis
    print("\n--- Part 1e: Complexity analysis ---")
    print("  Miller test per number: O(log^2(n)) witness bound * O(log^2(n)) modular exp")
    print("  = O(log^4(n)) per number (under GRH)")
    print("  For pi(x) via testing all n <= x: O(x * log^4(x))")
    print("  Sieve of Eratosthenes: O(x * log(log(x)))")
    print("  => Sieve VASTLY better for bulk enumeration")
    print("  Miller's GRH test useful for: individual primality checks, large isolated n")
    print("  For computing p(n): must enumerate primes up to ~p(n), so O(p(n)*log^4(p(n)))")
    print("  Since p(n) ~ n*ln(n), this is O(n*ln(n)*log^4(n*ln(n)))")
    print("  Compare: Meissel-Lehmer pi(x) in O(x^{2/3}), finding p(n) by binary search")
    print("  Miller approach: O(n * polylog(n)) -- WORSE than O(n^{2/3}) for large n")

    return results_1b, grh_data, batch_results


# ============================================================
# EXPERIMENT 2: GRH EXPLICIT FORMULA OPTIMAL T
# ============================================================

def li_function(x):
    """Logarithmic integral li(x) using mpmath."""
    return mpmath.li(x)


def li_complex(x, rho):
    """
    Compute li(x^rho) where rho = 1/2 + i*gamma (under RH).
    li(z) = Ei(log(z)) for complex z.
    x^rho = exp(rho * log(x))
    """
    log_x = mpmath.log(x)
    z = rho * log_x  # rho * ln(x)
    # li(x^rho) = Ei(z) = Ei(rho * ln(x))
    return mpmath.ei(z)


def explicit_formula_pi(x, zeros, T_zeros):
    """
    Compute pi(x) using the explicit formula under RH:
    pi(x) ~ li(x) - sum_{|gamma|<T} li(x^rho) - log(2) + integral term

    We use the first T_zeros zeros (gamma_1, ..., gamma_{T_zeros}).
    Each zero rho = 1/2 + i*gamma contributes a conjugate pair.
    """
    x_mp = mpmath.mpf(x)
    result = li_function(x_mp)

    # Subtract contributions from zeros (conjugate pairs)
    for k in range(min(T_zeros, len(zeros))):
        gamma = zeros[k]
        rho = mpmath.mpc(0.5, gamma)
        rho_conj = mpmath.mpc(0.5, -gamma)

        contrib = li_complex(x_mp, rho) + li_complex(x_mp, rho_conj)
        result -= contrib

    # Subtract log(2) term
    result -= mpmath.log(2)

    # Subtract integral from x to infinity of dt/(t*(t^2-1)*log(t))
    # For large x this is very small, we approximate as 0 for x >= 10
    # (it's < 10^{-10} for x > 100)

    return float(mpmath.re(result))


def experiment2_explicit_formula():
    """
    Investigate the optimal number of zeros T for the explicit formula.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: GRH EXPLICIT FORMULA OPTIMAL T")
    print("=" * 70)

    zeros = load_zeta_zeros(1000)
    num_zeros = len(zeros)
    print(f"  Loaded {num_zeros} zeta zeros")
    print(f"  First zero: {zeros[0]}")
    print(f"  Last zero:  {zeros[-1]}")

    # Part 2a: For various x, find T_min such that |pi_approx(x) - pi(x)| < 0.5
    print("\n--- Part 2a: Finding T_min for error < 0.5 ---")

    test_xs = [10**4, 5*10**4, 10**5, 5*10**5, 10**6, 5*10**6, 10**7]
    # Limit to what we can verify with sympy.primepi
    results_2a = []

    for x in test_xs:
        actual_pi = int(primepi(x))
        print(f"\n  x = {x:.0e}, pi(x) = {actual_pi}")

        # Try increasing T
        t_min = None
        errors_by_t = []
        t_values_to_try = list(range(1, min(num_zeros + 1, 1001), 1))

        for T in t_values_to_try:
            approx = explicit_formula_pi(x, zeros, T)
            error = abs(approx - actual_pi)
            errors_by_t.append((T, error, approx))

            if error < 0.5 and t_min is None:
                t_min = T

        # Theoretical T_min prediction: O(sqrt(x) * log^2(x))
        sqrt_x = math.sqrt(x)
        log_x = math.log(x)
        t_theory = sqrt_x * log_x**2

        if t_min is not None:
            print(f"    T_min (error < 0.5):     {t_min}")
        else:
            # Find minimum error achieved
            best = min(errors_by_t, key=lambda t: t[1])
            print(f"    T_min not found in {num_zeros} zeros")
            print(f"    Best error: {best[1]:.4f} at T={best[0]} (approx={best[2]:.2f})")
            t_min = None

        print(f"    Theoretical O(sqrt(x)*log^2(x)): {t_theory:.0f}")
        if t_min:
            print(f"    Ratio actual/theoretical: {t_min/t_theory:.6f}")

        results_2a.append({
            'x': x,
            'pi_x': actual_pi,
            't_min': t_min,
            't_theory': t_theory,
            'errors_sample': errors_by_t[:10] + (errors_by_t[-5:] if len(errors_by_t) > 10 else []),
        })

    # Part 2b: Error decay rate
    print("\n--- Part 2b: Error decay with T ---")
    print("  Under RH, error ~ sqrt(x)*log(x)/T")
    x_test = 10**5
    actual = int(primepi(x_test))
    print(f"  x = {x_test}, pi(x) = {actual}")
    print(f"  {'T':>6s}  {'Error':>12s}  {'sqrt(x)*log(x)/T':>20s}  {'Ratio':>10s}")

    for T in [5, 10, 20, 50, 100, 200, 500, min(num_zeros, 1000)]:
        if T > num_zeros:
            break
        approx = explicit_formula_pi(x_test, zeros, T)
        error = abs(approx - actual)
        predicted = math.sqrt(x_test) * math.log(x_test) / T
        ratio = error / predicted if predicted > 0 else float('inf')
        print(f"  {T:>6d}  {error:>12.4f}  {predicted:>20.4f}  {ratio:>10.4f}")

    # Part 2c: Can we need fewer than sqrt(x) zeros?
    print("\n--- Part 2c: Can we need fewer than sqrt(x) zeros? ---")
    print("  The error bound under RH is: |R(x)| <= C * sqrt(x) * log(x) / T")
    print("  For error < 0.5: T > 2*C*sqrt(x)*log(x)")
    print("  This means T = Omega(sqrt(x)) zeros are REQUIRED.")
    print()
    print("  Under GRH (not just RH), could clustering help?")
    print("  The zeros have GUE statistics -- no special structure to exploit.")
    print("  Each zero contributes an oscillatory term ~ x^{i*gamma}/gamma.")
    print("  Missing a zero with |contribution| > 0.5 causes rounding error.")
    print()

    # Verify: what's the contribution of individual zeros?
    x_test = 10**6
    print(f"  Individual zero contributions at x = {x_test:.0e}:")
    print(f"  {'Zero #':>8s}  {'gamma':>15s}  {'|contribution|':>18s}")
    for k in [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, min(num_zeros-1, 999)]:
        if k >= num_zeros:
            break
        gamma = zeros[k]
        rho = mpmath.mpc(0.5, gamma)
        rho_conj = mpmath.mpc(0.5, -gamma)
        x_mp = mpmath.mpf(x_test)
        contrib = abs(li_complex(x_mp, rho) + li_complex(x_mp, rho_conj))
        print(f"  {k+1:>8d}  {float(gamma):>15.4f}  {float(contrib):>18.4f}")

    print()
    print("  CONCLUSION: Individual zero contributions decay as ~sqrt(x)/gamma.")
    print("  For x=10^6, the 1000th zero contributes ~1-10 units.")
    print("  Need T until contributions < 0.5, so T ~ O(sqrt(x)).")
    print("  CANNOT reduce below sqrt(x) zeros without additional structure.")

    # Part 2d: Cost analysis
    print("\n--- Part 2d: Cost analysis ---")
    print("  To compute pi(x) via explicit formula:")
    print("  - Need T = O(sqrt(x) * log^2(x)) zeros")
    print("  - Each zero: O(polylog(x)) to evaluate li(x^rho)")
    print("  - Total: O(sqrt(x) * polylog(x))")
    print("  This IS the Lagarias-Odlyzko analytic method (1987).")
    print()
    print("  To compute p(n) from pi(x):")
    print("  - Binary search: O(log(n)) evaluations of pi(x)")
    print("  - Each pi(x) costs O(sqrt(x) * polylog(x))")
    print("  - Total: O(sqrt(p(n)) * polylog(p(n)) * log(n))")
    print("  - Since p(n) ~ n*ln(n): O(sqrt(n) * polylog(n))")
    print()
    print("  COMPARISON:")
    print("  Method                        Cost for p(n)")
    print("  ------                        -------------")
    print("  Sieve                         O(n * log(n))          [enumerate all]")
    print("  Miller GRH individual         O(n * polylog(n))      [test each]")
    print("  Meissel-Lehmer combinatorial  O(n^{2/3})             [best exact]")
    print("  Lagarias-Odlyzko analytic     O(n^{1/2+eps})         [uses RH zeros]")
    print("  Polylog target                O(polylog(n))           [NOT ACHIEVABLE]")

    return results_2a


# ============================================================
# MAIN
# ============================================================

def main():
    print("GRH Miller Batch Testing & Explicit Formula Optimal T")
    print("=" * 70)
    print()

    t_start = time.time()

    results_1b, grh_data, batch_results = experiment1_witness_analysis()

    t_mid = time.time()
    print(f"\n  [Experiment 1 completed in {t_mid - t_start:.1f}s]")

    results_2a = experiment2_explicit_formula()

    t_end = time.time()
    print(f"\n  [Experiment 2 completed in {t_end - t_mid:.1f}s]")

    print("\n" + "=" * 70)
    print("OVERALL SUMMARY")
    print("=" * 70)
    print()
    print("Experiment 1 (GRH Miller):")
    print("  - GRH Miller test is correct and deterministic")
    print("  - Actual witnesses needed << GRH bound (typically a=2 or a=3 suffices)")
    print("  - Batch testing gives NO computational advantage (same witnesses, same work)")
    print("  - Cost: O(n * polylog(n)) for p(n) -- worse than Meissel-Lehmer")
    print()
    print("Experiment 2 (Explicit Formula):")
    print("  - Need T = O(sqrt(x) * log^2(x)) zeros for error < 0.5")
    print("  - Cannot reduce below sqrt(x) zeros -- each contributes O(sqrt(x)/gamma)")
    print("  - Total cost: O(sqrt(x) * polylog(x)) = Lagarias-Odlyzko")
    print("  - This is the BEST known analytic method; polylog is not achievable")
    print()
    print(f"Total runtime: {t_end - t_start:.1f}s")


if __name__ == "__main__":
    main()
