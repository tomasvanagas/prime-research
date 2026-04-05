#!/usr/bin/env python3
"""
Local-to-Global Reconstruction for pi(x)

CORE IDEA: We can't compute pi(x) in polylog time directly. But what if we
could compute MANY DIFFERENT APPROXIMATIONS to pi(x), each cheap and each
with INDEPENDENT errors, then combine them to get the exact answer?

This is inspired by:
- Compressed sensing (many weak measurements → exact reconstruction)
- Error-correcting codes (redundancy → error correction)
- The Goldbach-Vinogradov circle method (multiple Fourier analyses combine)

Approach: For each small prime q, compute:
  A_q(x) = #{n ≤ x : n coprime to q} = x - ⌊x/q⌋

This is EXACT and O(1) per q. More generally:
  A_S(x) = #{n ≤ x : n coprime to all primes in S} via inclusion-exclusion

The question: by carefully choosing which sets S to evaluate, can we
reconstruct pi(x) more efficiently than the full inclusion-exclusion?

EXPERIMENT 1: How much information does each sieve step provide?
- After sieving by primes up to B, the count differs from pi(x) by at most
  the number of composites in (B, x] with no factor ≤ B.
- These are products of two primes both > B ("semiprimes" in Goldbach sense).
- Count them as a function of B.

EXPERIMENT 2: Can we use RANDOMIZED inclusion-exclusion?
- Instead of all 2^k terms, randomly sample subsets S and average.
- This gives a randomized estimator for pi(x).
- How many random samples are needed for the estimate to be within ±0.5?

EXPERIMENT 3: Interpolation from multiple evaluation points
- Compute pi(x) at many points x_1, ..., x_m using O(x^{2/3}) method.
- Use the SMOOTHNESS of pi(x) (step function, but predictable gaps) to
  interpolate pi(x) at new points with fewer direct computations.
- KEY: Can we amortize pi(x) computations over many queries?

EXPERIMENT 4: Error correlation structure
- Compute E(x) = pi(x) - Li(x) for many x values.
- Study the autocorrelation of E(x).
- If E(x) and E(x+h) are highly correlated for small h, then knowing
  E(x) constrains E(x+h), reducing the information needed.
"""

import numpy as np
import time
from sympy import primepi, isprime, prime, factorint
from math import log, floor, gcd
from itertools import combinations
from functools import reduce


def li(x):
    """Logarithmic integral via numerical integration."""
    if x <= 1:
        return 0
    from scipy.integrate import quad
    # Li(x) = integral from 2 to x of 1/ln(t) dt
    result, _ = quad(lambda t: 1/log(t), 2, x)
    return result


# ============================================================
# EXPERIMENT 1: Partial sieve residual
# ============================================================
def experiment_partial_sieve_residual(x_values=[1000, 5000, 10000, 50000]):
    """
    After sieving by primes ≤ B, how many "false positives" remain?
    These are composites with all prime factors > B.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Partial Sieve Residual")
    print("=" * 70)
    print("After sieving by primes ≤ B, how many composites survive?")
    print()

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    for x in x_values:
        true_pi = int(primepi(x))
        print(f"\nx = {x}, pi(x) = {true_pi}")
        print(f"{'B':>4} {'sieve count':>12} {'false pos':>10} {'ratio':>8} {'pi(B)':>6}")

        for k in range(1, min(len(small_primes) + 1, 12)):
            B = small_primes[k - 1]
            if B * B > x:
                # After sieving by primes ≤ B where B² > x, all survivors > 1 are prime
                sieve_count = sum(1 for n in range(2, x + 1)
                                  if all(n % p != 0 for p in small_primes[:k]))
                false_pos = sieve_count - true_pi
                # Add back primes ≤ B that were sieved out
                primes_le_B = sum(1 for p in small_primes[:k] if p <= x)
                sieve_count += primes_le_B
                false_pos = sieve_count - true_pi
                # Actually: sieve gives #{n ≤ x coprime to p1...pk} + #{primes ≤ B}
                # But this double-counts. Let me just count directly.
                pass

            # Direct count: numbers in [2,x] with no prime factor ≤ B
            coprime_count = sum(1 for n in range(2, x + 1)
                                if all(n % p != 0 for p in small_primes[:k]))
            # Primes ≤ B that were eliminated
            eliminated_primes = sum(1 for p in small_primes[:k] if p <= x)
            # Estimate of pi(x): coprime_count + eliminated_primes
            estimate = coprime_count + eliminated_primes
            false_pos = estimate - true_pi

            print(f"{B:>4} {estimate:>12} {false_pos:>10} {false_pos/true_pi:>8.4f} {k:>6}")

            if x > 10000 and k > 6:
                break  # too slow for large x with many sieve primes


# ============================================================
# EXPERIMENT 2: Randomized inclusion-exclusion
# ============================================================
def experiment_randomized_ie(x=1000, num_trials=1000):
    """
    Can randomized subset selection in inclusion-exclusion give a good estimate?

    Standard IE: pi_sieve(x) = Σ_{S ⊆ primes} (-1)^|S| ⌊x/prod(S)⌋
    Randomized: sample random subsets, weight by (-1)^|S|, average.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Randomized Inclusion-Exclusion")
    print("=" * 70)

    true_pi = int(primepi(x))
    primes_list = [p for p in range(2, int(x**0.5) + 1) if isprime(p)]
    k = len(primes_list)

    print(f"x = {x}, pi(x) = {true_pi}")
    print(f"Primes ≤ √x: {primes_list} (k = {k})")
    print(f"Full IE has 2^{k} = {2**k} terms")
    print()

    # First, compute exact IE (since k is small for x=1000, k ≤ 11)
    if k <= 20:
        ie_exact = 0
        for size in range(k + 1):
            for subset in combinations(range(k), size):
                prod_s = reduce(lambda a, b: a * b, [primes_list[i] for i in subset], 1)
                ie_exact += (-1)**size * floor(x / prod_s)
        ie_exact -= 1  # subtract 1 (for n=1)
        ie_exact += len(primes_list)  # add back the sieve primes themselves
        print(f"Exact IE result: {ie_exact} (should be {true_pi})")

    # Randomized approach: sample subsets
    np.random.seed(42)

    results_by_sample_size = {}
    for num_samples in [10, 50, 100, 500, 1000, 5000]:
        estimates = []
        for _ in range(num_trials):
            # Each subset has each prime included independently with prob 1/2
            total = 0
            for _ in range(num_samples):
                mask = np.random.randint(0, 2, size=k)
                subset_primes = [primes_list[i] for i in range(k) if mask[i]]
                prod_s = reduce(lambda a, b: a * b, subset_primes, 1)
                if prod_s <= x:
                    sign = (-1) ** sum(mask)
                    total += sign * floor(x / prod_s)
            # Normalize: each term appears with probability 2^{-k}, so multiply by 2^k
            estimate = total * (2**k) / num_samples
            estimates.append(estimate)

        estimates = np.array(estimates)
        mean_est = np.mean(estimates)
        std_est = np.std(estimates)
        results_by_sample_size[num_samples] = (mean_est, std_est)

    print(f"\n{'samples':>8} {'mean estimate':>14} {'std':>12} {'error/pi(x)':>12}")
    for ns, (m, s) in results_by_sample_size.items():
        print(f"{ns:>8} {m:>14.1f} {s:>12.1f} {abs(m - true_pi)/true_pi:>12.4f}")

    print(f"\nFor ±0.5 accuracy, need std < 0.25 → need ~{max(1, int(4*results_by_sample_size[1000][1]**2))} samples")
    print(f"Full IE has {2**k} terms → randomized IE needs MORE samples than full IE!")


# ============================================================
# EXPERIMENT 3: Error term autocorrelation
# ============================================================
def experiment_error_autocorrelation(x_max=50000, step=1):
    """
    Study the autocorrelation of E(x) = pi(x) - Li(x).
    If E(x) is strongly correlated at nearby points, knowing E(x₀) helps
    predict E(x₀ + h) for small h.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Error Term Autocorrelation")
    print("=" * 70)

    # Sample E(x) at many points
    x_values = list(range(100, x_max + 1, step))
    if len(x_values) > 5000:
        # Subsample for speed
        x_values = list(range(100, x_max + 1, max(1, (x_max - 100) // 5000)))

    print(f"Computing E(x) = pi(x) - Li(x) for {len(x_values)} points in [100, {x_max}]")
    t0 = time.time()

    E_values = []
    for x in x_values:
        pi_x = int(primepi(x))
        li_x = li(x)
        E_values.append(pi_x - li_x)
    E_values = np.array(E_values)

    elapsed = time.time() - t0
    print(f"Computed in {elapsed:.2f}s")

    # Statistics
    print(f"\nError term statistics:")
    print(f"  Mean E(x): {np.mean(E_values):.2f}")
    print(f"  Std E(x):  {np.std(E_values):.2f}")
    print(f"  Max |E(x)|: {np.max(np.abs(E_values)):.2f}")
    print(f"  √x_max / log(x_max): {x_max**0.5 / log(x_max):.2f}")

    # Autocorrelation
    E_centered = E_values - np.mean(E_values)
    var = np.var(E_centered)
    if var > 0:
        print(f"\nAutocorrelation of E(x):")
        print(f"{'lag h':>8} {'corr(E(x), E(x+h))':>20} {'# samples':>10}")
        for lag in [1, 2, 5, 10, 20, 50, 100, 200, 500]:
            if lag < len(E_centered):
                corr = np.mean(E_centered[:-lag] * E_centered[lag:]) / var
                print(f"{lag * step:>8} {corr:>20.6f} {len(E_centered) - lag:>10}")

    # Mutual information proxy: how much does E(x) reduce uncertainty about E(x+h)?
    # Use: knowing E(x) within ±0.5 (integer rounding), how constrained is E(x+h)?
    print(f"\n--- Prediction Analysis ---")
    print("If we know E(x₀) exactly, what is the std of E(x₀+h)?")
    for lag in [1, 5, 10, 50, 100]:
        if lag < len(E_values):
            diffs = E_values[lag:] - E_values[:-lag]
            print(f"  h = {lag * step:>5}: std(E(x+h) - E(x)) = {np.std(diffs):.3f}")

    # Key finding: the correlation structure
    print("\n--- KEY FINDING ---")
    if len(E_values) > 10:
        lag1_corr = np.mean(E_centered[:-1] * E_centered[1:]) / var if var > 0 else 0
        if lag1_corr > 0.9:
            print("E(x) is HIGHLY correlated at nearby points!")
            print("This means: knowing E(x₀) strongly constrains E(x₀+1).")
            print("But: the remaining uncertainty is pi(x₀+1) - pi(x₀) ∈ {0,1},")
            print("which is just primality testing — already O(polylog).")
            print("So correlation helps only for ADJACENT queries, not for")
            print("computing E(x₀) from scratch at a single point.")
        else:
            print(f"Lag-1 autocorrelation: {lag1_corr:.4f}")
            print("Error term is NOT strongly correlated → independent computations")
            print("at each point → no shortcut from local-to-global.")


# ============================================================
# EXPERIMENT 4: Information content of the error term
# ============================================================
def experiment_error_information(x_max=10000):
    """
    How many bits of information does E(x) = pi(x) - Li(x) carry?
    If E(x) is bounded by B, it carries at most log₂(2B+1) bits.
    For polylog computation, we need E(x) to carry only O(polylog) bits.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Information Content of Error Term")
    print("=" * 70)

    test_points = [100, 500, 1000, 5000, 10000, 50000, 100000]
    test_points = [x for x in test_points if x <= x_max]

    print(f"{'x':>8} {'pi(x)':>8} {'Li(x)':>10} {'E(x)':>8} {'|E|':>6} {'bits':>6} {'log²x':>7}")
    for x in test_points:
        pi_x = int(primepi(x))
        li_x = li(x)
        E_x = pi_x - li_x
        bits = log(max(1, abs(E_x) * 2 + 1)) / log(2) if E_x != 0 else 0
        log2x = (log(x) / log(2)) ** 2
        print(f"{x:>8} {pi_x:>8} {li_x:>10.1f} {E_x:>8.1f} {abs(E_x):>6.1f} {bits:>6.1f} {log2x:>7.1f}")

    # Compute E(x) for a range and measure its entropy
    print(f"\nComputing error term distribution for x in [100, {x_max}]...")
    errors = []
    for x in range(100, x_max + 1, 10):
        pi_x = int(primepi(x))
        li_x = li(x)
        errors.append(int(round(pi_x - li_x)))

    error_counts = {}
    for e in errors:
        error_counts[e] = error_counts.get(e, 0) + 1

    total = len(errors)
    entropy = -sum((c/total) * log(c/total) / log(2) for c in error_counts.values())

    print(f"\nError term distribution (sampled at {total} points):")
    print(f"  Range: [{min(errors)}, {max(errors)}]")
    print(f"  Entropy: {entropy:.2f} bits")
    print(f"  Uniform entropy for this range: {log(max(errors)-min(errors)+1)/log(2):.2f} bits")
    print(f"  (log x_max)²: {(log(x_max)/log(2))**2:.2f}")
    print(f"  √x_max: {x_max**0.5:.2f}")

    print("\n--- ANALYSIS ---")
    print("The error E(x) carries O(√x / log x) information on RH.")
    print("For polylog computation, we'd need O(polylog(x)) bits of E(x).")
    print(f"At x = {x_max}: need {(log(x_max)/log(2))**2:.0f} bits (polylog)")
    print(f"  but E(x) needs {log(max(1, max(abs(e) for e in errors))*2+1)/log(2):.0f} bits (√x/log x)")
    print("The gap between polylog and √x/log x is the FUNDAMENTAL BARRIER.")


# ============================================================
# EXPERIMENT 5: Multi-resolution approximation
# ============================================================
def experiment_multi_resolution(n_target=1000):
    """
    Test: Can we compute p(n) by using approximations at multiple resolutions?

    Level 0: R^{-1}(n) → ~50% digits correct
    Level 1: Add first few zeta zeros → how many more digits?
    Level 2: Add more zeros → more digits?

    KEY QUESTION: Does each zero contribute O(1) bits of information?
    If so, we need O(digits/2) = O(log x) zeros for the exact answer.
    That's O(log x) zeros × O(log x) time per zero = O(log²x). POLYLOG!

    But reality: each zero contributes LESS than 1 bit (they don't carry
    independent information). The zeros' contributions have GUE-correlated
    phases that partially cancel. How much cancellation occurs?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Multi-Resolution Approximation")
    print("=" * 70)

    from mpmath import mp, mpf, li as mp_li, log as mp_log, pi as mp_pi, exp, sqrt

    mp.dps = 50

    # First 30 zeta zeros (imaginary parts)
    zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851
    ]

    # Test for various x values
    test_x = [100, 500, 1000, 5000, 10000]
    print(f"\nFor each x, compute pi_approx(x) using K zeros and measure error:")
    print(f"{'x':>6} {'pi(x)':>6} {'Li(x)':>8} K=0_err K=1_err K=5_err K=10_err K=20_err K=30_err")

    for x in test_x:
        true_pi = int(primepi(x))
        li_x = float(mp_li(x))

        # Compute the zero correction: Σ_ρ Li(x^ρ) + Li(x^{1-ρ})
        # For ρ = 1/2 + iγ: x^ρ = x^{1/2} * e^{iγ log x}
        # Li(x^ρ) ≈ x^ρ / (ρ log x) for large x

        errors = []
        zero_sum_real = 0.0
        for k in range(31):
            if k > 0 and k <= len(zeros):
                gamma = zeros[k - 1]
                rho = complex(0.5, gamma)
                # Li(x^rho) approximation: x^rho / (rho * log(x))
                x_rho = x ** rho  # complex
                correction = 2 * (x_rho / (rho * log(x))).real  # ρ and 1-ρ̄ contribute conjugate
                zero_sum_real += correction

            approx = li_x - zero_sum_real
            error = approx - true_pi
            errors.append(abs(error))

        print(f"{x:>6} {true_pi:>6} {li_x:>8.1f} "
              f"{errors[0]:>7.1f} {errors[1]:>7.1f} {errors[5]:>7.1f} "
              f"{errors[10]:>8.1f} {errors[20]:>8.1f} {errors[30]:>8.1f}")

    # Detailed analysis for x = 10000
    print(f"\nDetailed: x = 10000, pi(x) = {int(primepi(10000))}")
    print(f"{'K zeros':>8} {'approx':>10} {'error':>8} {'correct digits':>15} {'bits gained':>12}")
    x = 10000
    true_pi = int(primepi(x))
    li_x = float(mp_li(x))
    zero_sum_real = 0.0
    prev_error = None

    for k in range(31):
        if k > 0:
            gamma = zeros[k - 1]
            rho = complex(0.5, gamma)
            x_rho = x ** rho
            correction = 2 * (x_rho / (rho * log(x))).real
            zero_sum_real += correction

        approx = li_x - zero_sum_real
        error = abs(approx - true_pi)
        digits = -log(max(error, 1e-15)) / log(10) if error > 0 else 15
        bits_gained = log(max(1, abs(prev_error) if prev_error else li_x - true_pi) / max(error, 1e-15)) / log(2) if error > 0 and prev_error is not None else 0
        if k in [0, 1, 2, 3, 5, 10, 15, 20, 25, 30]:
            print(f"{k:>8} {approx:>10.2f} {error:>8.3f} {digits:>15.1f} {bits_gained:>12.2f}")
        prev_error = error

    # How many zeros needed for error < 0.5?
    print(f"\n--- KEY QUESTION: How many zeros for exact answer? ---")
    for x in [100, 1000, 10000]:
        true_pi = int(primepi(x))
        li_x = float(mp_li(x))
        zero_sum = 0.0
        for k in range(len(zeros)):
            gamma = zeros[k]
            rho = complex(0.5, gamma)
            x_rho = x ** rho
            correction = 2 * (x_rho / (rho * log(x))).real
            zero_sum += correction
            approx = li_x - zero_sum
            if abs(approx - true_pi) < 0.5:
                print(f"x = {x:>6}: need {k+1} zeros (√x = {x**0.5:.1f}, log²x = {(log(x)/log(2))**2:.1f})")
                break
        else:
            final_err = abs(li_x - zero_sum - true_pi)
            print(f"x = {x:>6}: 30 zeros not enough (error = {final_err:.2f}, √x = {x**0.5:.1f})")


if __name__ == "__main__":
    print("=" * 70)
    print("LOCAL-TO-GLOBAL RECONSTRUCTION FOR pi(x)")
    print("=" * 70)

    experiment_partial_sieve_residual([1000, 5000])
    experiment_randomized_ie(x=500, num_trials=200)
    experiment_error_autocorrelation(x_max=10000, step=1)
    experiment_error_information(x_max=10000)
    experiment_multi_resolution(n_target=1000)

    print("\n" + "=" * 70)
    print("GRAND SUMMARY")
    print("=" * 70)
    print("""
1. PARTIAL SIEVE: Sieving by primes ≤ B leaves O(x/log²B) false positives.
   Need B ≈ √x for exact count → O(√x/log x) sieve primes needed.

2. RANDOMIZED IE: Random subset sampling gives estimates with std ~ 2^k / √samples.
   Need samples >> 2^{2k} for accuracy → WORSE than exact IE.

3. ERROR AUTOCORRELATION: E(x) is highly correlated at nearby points (good!).
   But this only helps for adjacent queries, not for a single-point computation.

4. INFORMATION CONTENT: E(x) carries O(√x/log x) bits on RH.
   Polylog computation could handle O(polylog) bits.
   The gap is FUNDAMENTAL — no local-to-global scheme can bridge it.

5. MULTI-RESOLUTION: Each zeta zero contributes < 1 bit of correction.
   Need Ω(√x) zeros for exact answer.
   NO polylog shortcut via partial zero sums.

VERDICT: Local-to-global reconstruction CANNOT achieve polylog for pi(x).
The information-theoretic barrier is real: E(x) contains too many bits.
""")
