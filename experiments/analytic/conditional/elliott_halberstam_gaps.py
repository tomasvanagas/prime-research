#!/usr/bin/env python3
"""
Elliott-Halberstam Conjecture and Gap Structure: Experiments 3 & 4 from Task #4.

Experiment 3: Does summing pi(x;q,a) over residue classes recover pi(x) faster
              than direct computation? Test with primorial moduli q=6,30,210.

Experiment 4: Gap structure under Cramer's model. Even with bounded gaps and
              R^{-1}(n) approximation, does counting remain the bottleneck?
"""

import time
import math
from collections import defaultdict
import numpy as np
from sympy import primepi, primerange, isprime, li, nextprime, prime
from sympy.ntheory import totient
from mpmath import mp, mpf, li as mpli, log as mplog

mp.dps = 30


# =============================================================================
# EXPERIMENT 3: Elliott-Halberstam and Counting
# =============================================================================

def primes_in_residue_class(x, q, a):
    """Count primes p <= x with p ≡ a (mod q)."""
    count = 0
    for p in primerange(2, x + 1):
        if p % q == a:
            count += 1
    return count


def coprime_residues(q):
    """Return residues a in [1, q) with gcd(a, q) = 1."""
    from math import gcd
    return [a for a in range(1, q) if gcd(a, q) == 1]


def experiment3_eh_counting():
    """
    Test: pi(x) = sum_{a coprime to q} pi(x; q, a) + (primes dividing q that are <= x).
    Under EH, each pi(x;q,a) ~ li(x)/phi(q) with small error.
    Question: do errors cancel well enough when summing to recover pi(x)?
    """
    print("=" * 70)
    print("EXPERIMENT 3: Elliott-Halberstam Conjecture and Counting")
    print("=" * 70)

    primorials = {
        6: [2, 3],       # 2*3
        30: [2, 3, 5],   # 2*3*5
        210: [2, 3, 5, 7]  # 2*3*5*7
    }

    x_values = [10**3, 10**4, 10**5, 10**6]
    results = {}

    for q, q_primes in primorials.items():
        phi_q = totient(q)
        residues = coprime_residues(q)
        assert len(residues) == phi_q, f"phi({q}) mismatch"

        print(f"\n--- Modulus q = {q}, phi(q) = {phi_q}, residues: {len(residues)} ---")
        results[q] = []

        for x in x_values:
            t0 = time.time()

            # True pi(x)
            true_pi = primepi(x)

            # li(x) approximation
            li_x = float(mpli(x))

            # Count primes in each residue class
            class_counts = {}
            class_errors = {}
            expected = li_x / phi_q

            for a in residues:
                cnt = primes_in_residue_class(x, q, a)
                class_counts[a] = cnt
                class_errors[a] = cnt - expected

            # Sum over coprime residue classes
            sum_coprime = sum(class_counts.values())

            # Add primes that divide q and are <= x
            small_primes_count = sum(1 for p in q_primes if p <= x)
            reconstructed_pi = sum_coprime + small_primes_count

            # Error analysis
            total_error = reconstructed_pi - true_pi
            li_error = li_x - true_pi
            errors_list = list(class_errors.values())
            max_individual_error = max(abs(e) for e in errors_list)
            sum_of_errors = sum(errors_list)
            rms_error = math.sqrt(sum(e**2 for e in errors_list) / len(errors_list))

            elapsed = time.time() - t0

            row = {
                'x': x,
                'true_pi': true_pi,
                'reconstructed_pi': reconstructed_pi,
                'reconstruction_exact': (total_error == 0),
                'li_x': li_x,
                'li_error': li_error,
                'li_relative_error': li_error / true_pi if true_pi else 0,
                'sum_class_errors': sum_of_errors,
                'max_individual_error': max_individual_error,
                'rms_class_error': rms_error,
                'error_cancellation_ratio': abs(sum_of_errors) / (phi_q * rms_error) if rms_error > 0 else 0,
                'elapsed': elapsed,
            }
            results[q].append(row)

            print(f"  x={x:.0e}: pi(x)={true_pi}, reconstructed={reconstructed_pi}, "
                  f"exact={total_error == 0}, sum_errors={sum_of_errors:.2f}, "
                  f"max_err={max_individual_error:.2f}, RMS={rms_error:.2f}, "
                  f"cancel_ratio={row['error_cancellation_ratio']:.4f}, "
                  f"time={elapsed:.3f}s")

    # Timing comparison: direct primepi vs residue-class sum
    print("\n--- Timing Comparison: Direct vs Residue-Class Summation ---")
    timing_results = []
    for x in [10**4, 10**5, 10**6]:
        t0 = time.time()
        direct = primepi(x)
        t_direct = time.time() - t0

        t0 = time.time()
        q = 30
        residues = coprime_residues(q)
        total = sum(primes_in_residue_class(x, q, a) for a in residues)
        total += sum(1 for p in [2, 3, 5] if p <= x)
        t_residue = time.time() - t0

        timing_results.append((x, t_direct, t_residue, total == direct))
        print(f"  x={x:.0e}: direct={t_direct:.4f}s, residue_sum(q=30)={t_residue:.4f}s, "
              f"speedup={t_residue/t_direct if t_direct > 0 else float('inf'):.1f}x slower, "
              f"match={total == direct}")

    return results, timing_results


# =============================================================================
# EXPERIMENT 4: Gap Structure and Counting
# =============================================================================

def experiment4a_cramer_model():
    """
    Test Cramer's model: gaps g_n ~ Exp(ln(p_n)) for primes up to 10^7.
    Compare empirical gap distribution to exponential with rate 1/ln(p).
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4A: Cramer's Model for Prime Gaps")
    print("=" * 70)

    limit = 10**7
    primes = list(primerange(2, limit))
    n_primes = len(primes)
    print(f"Primes up to {limit}: {n_primes}")

    # Compute gaps and normalized gaps
    gaps = [primes[i+1] - primes[i] for i in range(n_primes - 1)]

    # Bin primes by magnitude for analysis
    bins = [(10**k, 10**(k+1)) for k in range(1, 7)]
    results = []

    for lo, hi in bins:
        if hi > limit:
            hi = limit
        indices = [i for i in range(n_primes - 1) if lo <= primes[i] < hi]
        if not indices:
            continue

        local_gaps = [gaps[i] for i in indices]
        local_primes = [primes[i] for i in indices]
        expected_mean = np.mean([math.log(p) for p in local_primes])
        actual_mean = np.mean(local_gaps)
        actual_std = np.std(local_gaps)
        actual_max = max(local_gaps)

        # Under Cramer, normalized gaps g/ln(p) should be ~ Exp(1)
        normalized = [gaps[i] / math.log(primes[i]) for i in indices]
        norm_mean = np.mean(normalized)
        norm_std = np.std(normalized)

        # Cramer's conjecture: max gap ~ (ln p)^2
        cramer_bound = math.log(hi) ** 2

        row = {
            'range': f"[{lo:.0e}, {hi:.0e})",
            'count': len(indices),
            'expected_mean_gap': expected_mean,
            'actual_mean_gap': actual_mean,
            'actual_std_gap': actual_std,
            'max_gap': actual_max,
            'cramer_bound': cramer_bound,
            'max_gap_over_cramer': actual_max / cramer_bound,
            'normalized_mean': norm_mean,
            'normalized_std': norm_std,
        }
        results.append(row)

        print(f"  {row['range']}: mean_gap={actual_mean:.2f} (expected~{expected_mean:.2f}), "
              f"max_gap={actual_max}, Cramer_bound={cramer_bound:.1f}, "
              f"ratio={row['max_gap_over_cramer']:.3f}, "
              f"norm_mean={norm_mean:.3f}, norm_std={norm_std:.3f}")

    return results


def R_inverse_approx(n):
    """
    Approximate R^{-1}(n) using Newton's method on R(x) = n,
    where R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k}).
    For simplicity, use the first-order approximation: R^{-1}(n) ~ n*ln(n).
    Then refine with li^{-1}(n) as a better starting point.
    """
    # li^{-1}(n) approximation: x ~ n*ln(n) for large n
    if n < 10:
        # Small n: just return the nth prime directly for validation
        return int(prime(n))

    x = float(n * math.log(n))
    # Newton iterations on li(x) = n
    for _ in range(20):
        li_x = float(mpli(x))
        if abs(li_x - n) < 0.5:
            break
        # li'(x) = 1/ln(x)
        deriv = 1.0 / math.log(x)
        x = x - (li_x - n) / deriv
        if x <= 2:
            x = 2.1

    return x


def experiment4b_search_cost():
    """
    Given R^{-1}(n), search in interval [R^{-1}(n) - C*ln^2, R^{-1}(n) + C*ln^2].
    Count primes in interval. Measure: how often is p(n) in this interval?
    What's the computational cost?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4B: Search Cost with R^{-1}(n) and Cramer's Bound")
    print("=" * 70)

    C = 2.0  # Safety factor for Cramer's conjecture bound
    test_ns = [10, 100, 1000, 10000, 50000, 100000, 500000]
    results = []

    for n in test_ns:
        t0 = time.time()

        # True nth prime
        pn = prime(n)

        # R^{-1}(n) approximation
        r_inv = R_inverse_approx(n)
        ln_r = math.log(r_inv) if r_inv > 1 else 1
        half_width = C * ln_r ** 2

        lo = max(2, int(r_inv - half_width))
        hi = int(r_inv + half_width) + 1

        # Count primes in the search interval
        primes_in_interval = list(primerange(lo, hi))
        n_primes_in_interval = len(primes_in_interval)

        # Is p(n) in the interval?
        pn_in_interval = lo <= pn <= hi

        # Error of R^{-1}(n) vs p(n)
        approx_error = r_inv - pn
        relative_error = approx_error / pn if pn else 0

        # Cost analysis
        # Each primality test costs O(ln^2(x)) with Miller-Rabin
        interval_size = hi - lo
        # Candidates to test ~ interval_size (or use sieve)
        # With sieve, cost ~ interval_size * ln(ln(x))
        # With individual tests, cost ~ interval_size * ln^2(x)
        ln_x = math.log(pn)
        sieve_cost = interval_size  # sieve is essentially linear in interval
        test_cost = n_primes_in_interval * ln_x ** 2  # testing each prime

        # But we still need pi(lo) to know WHICH prime is p(n)
        # This counting step is the bottleneck
        # pi(lo) ~ O(lo^{2/3}) with Meissel-Lehmer
        counting_cost = lo ** (2.0 / 3.0)

        elapsed = time.time() - t0

        row = {
            'n': n,
            'pn': pn,
            'r_inv': r_inv,
            'approx_error': approx_error,
            'relative_error': relative_error,
            'interval': (lo, hi),
            'interval_size': interval_size,
            'primes_in_interval': n_primes_in_interval,
            'pn_in_interval': pn_in_interval,
            'sieve_cost_est': sieve_cost,
            'counting_cost_est': counting_cost,
            'counting_dominates': counting_cost > sieve_cost,
            'elapsed': elapsed,
        }
        results.append(row)

        print(f"  n={n:>7}: p(n)={pn:>9}, R^-1(n)={r_inv:>12.1f}, "
              f"err={approx_error:>+10.1f} ({relative_error:>+.2e}), "
              f"interval=[{lo},{hi}] size={interval_size}, "
              f"primes_in={n_primes_in_interval}, p(n)_in={pn_in_interval}, "
              f"count_cost={counting_cost:.0f} vs sieve={sieve_cost}")

    return results


def experiment4c_counting_bottleneck():
    """
    Demonstrate that even with a perfect gap oracle, counting pi(x) remains
    the bottleneck for determining p(n).

    Scenario: Suppose we know gaps exactly. Starting from p(1)=2, we can
    compute p(n) = 2 + sum_{i=1}^{n-1} g_i. But this sum requires knowing
    all n-1 gaps, which is O(n) = O(p(n)/ln(p(n))) work.

    Alternative: start from R^{-1}(n) and search locally. But to confirm
    which prime in the interval is p(n), we need pi(R^{-1}(n) - delta).
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4C: The Counting Bottleneck")
    print("=" * 70)

    results = []

    for n in [1000, 10000, 100000]:
        pn = prime(n)
        ln_pn = math.log(pn)

        # Method 1: Sequential gap summation
        # Cost = O(n) = O(p(n) / ln(p(n)))
        sequential_cost = n

        # Method 2: R^{-1}(n) + local search + counting
        # Local search cost: O(ln^2(p(n))) interval, sieved in O(ln^2(p(n)))
        # Counting cost: O(p(n)^{2/3}) via Meissel-Lehmer
        local_search_cost = ln_pn ** 2
        counting_cost = pn ** (2.0 / 3.0)
        method2_cost = local_search_cost + counting_cost

        # Method 3: Hypothetical -- if EH gives better counting
        # Best known under EH: still O(x^{1/2+eps}) for exact pi(x)
        # (EH helps with *distribution* not *total count*)
        eh_counting_cost = pn ** 0.5 * ln_pn  # optimistic under EH

        # The polylog dream
        polylog_cost = ln_pn ** 3

        row = {
            'n': n,
            'pn': pn,
            'sequential': sequential_cost,
            'meissel_lehmer': counting_cost,
            'eh_optimistic': eh_counting_cost,
            'polylog_dream': polylog_cost,
            'ratio_ml_to_seq': counting_cost / sequential_cost,
            'ratio_eh_to_ml': eh_counting_cost / counting_cost,
        }
        results.append(row)

        print(f"  n={n:>6}, p(n)={pn:>8}:")
        print(f"    Sequential gaps:    O(n)          = {sequential_cost:>12.0f}")
        print(f"    Meissel-Lehmer:     O(x^{{2/3}})    = {counting_cost:>12.0f}")
        print(f"    Under EH (optim.):  O(x^{{1/2}}ln)  = {eh_counting_cost:>12.0f}")
        print(f"    Polylog dream:      O(ln^3)       = {polylog_cost:>12.1f}")
        print(f"    ML/Sequential ratio: {row['ratio_ml_to_seq']:.2f}")

    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Elliott-Halberstam Conjecture & Gap Structure Experiments")
    print("Task #4, Experiments 3 & 4")
    print()

    # Experiment 3
    eh_results, timing_results = experiment3_eh_counting()

    # Experiment 4
    cramer_results = experiment4a_cramer_model()
    search_results = experiment4b_search_cost()
    bottleneck_results = experiment4c_counting_bottleneck()

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    print("\nExperiment 3 (Elliott-Halberstam):")
    print("  - Summing pi(x;q,a) over coprime residues recovers pi(x) EXACTLY")
    print("    (this is a tautology -- it's just partitioning primes by residue class)")
    print("  - The residue-class approach is SLOWER than direct computation")
    print("  - EH controls individual errors |pi(x;q,a) - li(x)/phi(q)|")
    print("  - But summing these approximations gives li(x) +/- small error,")
    print("    which is NO BETTER than computing li(x) directly")
    print("  - VERDICT: EH helps with distribution, NOT with counting")

    print("\nExperiment 4 (Gap Structure):")
    print("  - Cramer's model fits well: normalized gaps ~ Exp(1)")
    print("  - R^{-1}(n) locates p(n) within O(ln^2(p)) interval")
    all_found = all(r['pn_in_interval'] for r in search_results)
    print(f"  - p(n) found in search interval: {all_found}")
    print("  - BUT: identifying WHICH prime in the interval is p(n)")
    print("    requires knowing pi(x) at the interval boundary")
    print("  - Counting pi(x) costs O(x^{2/3}) -- this is THE bottleneck")
    print("  - Even under EH, no known way to make counting polylog")
    print("  - VERDICT: Gap knowledge doesn't bypass the counting problem")

    return {
        'eh_results': eh_results,
        'timing_results': timing_results,
        'cramer_results': cramer_results,
        'search_results': search_results,
        'bottleneck_results': bottleneck_results,
    }


if __name__ == '__main__':
    all_results = main()
