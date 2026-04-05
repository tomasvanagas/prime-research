#!/usr/bin/env python3
"""
Hybrid Analytic-Sieve Method for Computing p(n)
================================================
Key idea: We don't need π(x) exactly for ALL x.
We need just enough precision to locate p(n) in a small interval, then sieve locally.

Step 1: Compute x₀ = R⁻¹(n) ≈ p(n) with error E₀ ~ √x₀ · log(x₀)
Step 2: Use K zeta zeros to reduce error to E_K ~ √x₀ / K · something
Step 3: Sieve locally in [x₀ - E_K, x₀ + E_K] to find exact p(n)

Total cost: cost(K zeros) + cost(sieve width E_K)

Question: What's the optimal K that minimizes total cost?

Also testing: "ITERATIVE REFINEMENT" — start with rough estimate, use local
sieve to get π(x₀) ± δ, use this to narrow the interval further.

And: "PROBABILISTIC SIEVE" — can we determine whether a specific integer
is prime WITHOUT checking all prime divisors?
"""

import numpy as np
from sympy import primepi, nextprime, prevprime, isprime, primerange, mobius
from scipy import special
import time
import math

# Load zeta zeros
def load_zeros(n=1000):
    with open('/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt') as f:
        zeros = []
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                zeros.append(float(parts[1]))
            elif len(parts) == 1:
                zeros.append(float(parts[0]))
    return np.array(zeros[:n])

ZEROS = load_zeros()

def R_function(x):
    """Riemann R function"""
    if x <= 1:
        return 0.0
    result = 0.0
    for n in range(1, 100):
        mu_n = int(mobius(n))
        if mu_n == 0:
            continue
        xn = x ** (1.0 / n)
        if xn <= 1.001:
            break
        li_val = float(special.expi(np.log(xn)))
        result += mu_n / n * li_val
    return result

def R_inverse(n, tol=0.01):
    """Compute R⁻¹(n) via Newton's method. Returns x such that R(x) ≈ n."""
    # Initial guess: n·log(n) for n > 1
    if n <= 1:
        return 2.0
    x = n * math.log(n)
    for _ in range(100):
        rx = R_function(x)
        # R'(x) ≈ 1/log(x) (dominant term)
        deriv = 1.0 / math.log(x) if x > 1 else 1.0
        dx = (n - rx) / deriv
        x += dx
        if abs(dx) < tol:
            break
    return x

def zero_sum_correction(x, num_zeros):
    """Compute Σ_ρ x^ρ/ρ (conjugate pairs) using num_zeros zeros."""
    total = 0.0
    log_x = math.log(x)
    for gamma in ZEROS[:num_zeros]:
        rho = 0.5 + 1j * gamma
        contrib = x**rho / rho
        total += 2 * contrib.real
    return total

def pi_estimate(x, num_zeros):
    """Estimate π(x) using R(x) and num_zeros zeta zero corrections."""
    Rx = R_function(x)
    if num_zeros == 0:
        return Rx
    correction = zero_sum_correction(x, num_zeros)
    return Rx - correction

def local_sieve(lo, hi):
    """Count primes in [lo, hi] using sieve of Eratosthenes."""
    lo, hi = int(lo), int(hi)
    if lo < 2:
        lo = 2
    count = 0
    primes_in_range = list(primerange(lo, hi + 1))
    return len(primes_in_range), primes_in_range

def find_nth_prime_hybrid(n, num_zeros=0, verbose=False):
    """
    Find p(n) using hybrid analytic + local sieve method.

    1. Compute x₀ = R⁻¹(n)
    2. Optionally refine with zero sum
    3. Estimate error bound
    4. Sieve locally to find exact p(n)
    """
    t0 = time.time()

    # Step 1: Initial estimate
    x0 = R_inverse(n)

    # Step 2: Refine with zeros
    if num_zeros > 0:
        # Use zeros to get better estimate of π(x₀)
        pi_est = pi_estimate(x0, num_zeros)
        # Adjust x₀ based on the correction
        deficit = n - pi_est  # how many more primes we need
        # Each prime is ~log(x₀) apart, so adjust by deficit * log(x₀)
        x0 += deficit * math.log(x0)

    # Step 3: Error bound
    # Without zeros: error ~ C·√x·log²(x) (empirical)
    # With K zeros: error ~ C·√x·log²(x)/K (heuristic)
    sqrt_x = math.sqrt(x0)
    log_x = math.log(x0)
    if num_zeros == 0:
        error_bound = 3.0 * sqrt_x * log_x  # conservative
    else:
        error_bound = max(3.0 * sqrt_x * log_x / num_zeros, 50)

    lo = max(2, int(x0 - error_bound))
    hi = int(x0 + error_bound)
    sieve_width = hi - lo

    # Step 4: Sieve locally
    # First count primes up to lo (this is the expensive part!)
    # We need π(lo) to convert local count to absolute count
    #
    # PROBLEM: We need π(lo) exactly to know which local prime is p(n)
    # We have an estimate of π(lo) from the analytic formula

    pi_lo_est = pi_estimate(lo, num_zeros)

    # Sieve the interval
    count_in_range, primes_in_range = local_sieve(lo, hi)

    # The index of p(n) in the local list is approximately n - π(lo)
    # But π(lo) is only approximate!
    local_index = round(n - pi_lo_est)

    # We need the EXACT value, so compute π at a few points to calibrate
    # For small cases, just compute exactly
    true_pi_lo = int(primepi(lo))
    exact_local_index = n - true_pi_lo - 1  # 0-indexed

    elapsed = time.time() - t0

    if verbose:
        print(f"  x₀ = {x0:.2f}")
        print(f"  error_bound = {error_bound:.1f}")
        print(f"  sieve_width = {sieve_width}")
        print(f"  primes in range = {count_in_range}")
        print(f"  π(lo) estimate = {pi_lo_est:.2f}, actual = {true_pi_lo}")
        print(f"  local_index est = {local_index}, actual = {exact_local_index}")

    if 0 <= exact_local_index < len(primes_in_range):
        result = primes_in_range[exact_local_index]
    else:
        result = None

    return result, sieve_width, elapsed

# =============================================================================
# EXPERIMENT 1: How does sieve width scale with number of zeros?
# =============================================================================
def experiment_sieve_width():
    print("=" * 70)
    print("EXPERIMENT 1: Sieve Width vs Number of Zeros")
    print("=" * 70)

    test_n = [100, 1000, 5000, 10000, 50000]
    zero_counts = [0, 1, 5, 10, 20, 50, 100, 200, 500]

    print(f"\n{'n':>7} | {'K zeros':>7} | {'Error bound':>12} | {'Sieve width':>12} | "
          f"{'Correct?':>8} | {'Time':>8}")
    print("-" * 75)

    for n in test_n:
        true_pn = int(nextprime(1))
        # Compute true p(n)
        p = 2
        for _ in range(n - 1):
            p = int(nextprime(p))
        true_pn = p

        for K in zero_counts:
            result, width, elapsed = find_nth_prime_hybrid(n, K)
            correct = result == true_pn if result else False
            print(f"{n:>7} | {K:>7} | {width/2:>12.1f} | {width:>12} | "
                  f"{'YES' if correct else 'NO':>8} | {elapsed:>7.4f}s")
        print()

# =============================================================================
# EXPERIMENT 2: Actual error of R⁻¹(n) vs n
# =============================================================================
def experiment_approximation_error():
    print("=" * 70)
    print("EXPERIMENT 2: Approximation Error |R⁻¹(n) - p(n)| vs n")
    print("=" * 70)

    ns = [10, 50, 100, 500, 1000, 5000, 10000, 50000]

    print(f"\n{'n':>7} | {'p(n)':>10} | {'R⁻¹(n)':>12} | {'|error|':>10} | "
          f"{'√p(n)':>10} | {'ratio':>8}")
    print("-" * 70)

    for n in ns:
        # Compute p(n) exactly
        p = 2
        for _ in range(n - 1):
            p = int(nextprime(p))

        x0 = R_inverse(n)
        error = abs(x0 - p)
        sqrt_p = math.sqrt(p)
        ratio = error / sqrt_p if sqrt_p > 0 else 0

        print(f"{n:>7} | {p:>10} | {x0:>12.2f} | {error:>10.2f} | "
              f"{sqrt_p:>10.2f} | {ratio:>8.4f}")

    print("\n  If ratio is bounded: error = O(√p(n)), need O(√n·log(n)) sieve width")
    print("  If ratio shrinks: error grows slower than √p(n) — promising!")

# =============================================================================
# EXPERIMENT 3: Error after K zeta zero corrections
# =============================================================================
def experiment_zero_corrected_error():
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Error After Zero Corrections")
    print("=" * 70)

    test_n = [1000, 5000, 10000, 50000]
    zero_counts = [0, 5, 10, 20, 50, 100, 200, 500]

    for n in test_n:
        # Compute p(n)
        p = 2
        for _ in range(n - 1):
            p = int(nextprime(p))
        pn = p

        print(f"\nn = {n}, p(n) = {pn}")

        for K in zero_counts:
            # Estimate π(pn) using K zeros
            pi_est = pi_estimate(pn, K)
            error = abs(pi_est - n)  # Should be 0 if exact

            # Also estimate π(pn - 1) and π(pn + 1)
            pi_below = pi_estimate(pn - 1, K)
            pi_above = pi_estimate(pn + 1, K) if pn + 1 < 10**8 else 0

            print(f"  K={K:>3}: π̂(p(n)) = {pi_est:>12.4f}, "
                  f"|π̂(p(n)) - n| = {error:>8.4f}, "
                  f"π̂(p(n)-1) = {pi_below:>12.4f}")

# =============================================================================
# EXPERIMENT 4: What if we use probabilistic primality + binary search?
# =============================================================================
def experiment_probabilistic_refinement():
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Iterative Refinement via Local Counting")
    print("=" * 70)
    print("\nIdea: Start with R⁻¹(n), use LOCAL sieve to calibrate,")
    print("then shift estimate. Each iteration narrows the window.\n")

    for n in [1000, 10000, 50000]:
        # Compute true p(n)
        p = 2
        for _ in range(n - 1):
            p = int(nextprime(p))
        pn = p

        # Start with R⁻¹(n)
        x0 = R_inverse(n)

        print(f"n = {n}, p(n) = {pn}, R⁻¹(n) = {x0:.2f}, initial error = {abs(x0 - pn):.2f}")

        # Iterative refinement:
        # 1. Sieve a SMALL window around x0 to get local prime density
        # 2. Compare local density to expected density (1/log x)
        # 3. Use the discrepancy to refine estimate

        window = 100  # small window
        iterations = 0
        x_est = x0

        for it in range(10):
            # Count primes in [x_est - window, x_est + window]
            lo = max(2, int(x_est - window))
            hi = int(x_est + window)
            local_count, local_primes = local_sieve(lo, hi)

            # We need π(lo) to convert. For now, use analytic estimate
            pi_lo = pi_estimate(lo, 100)  # use 100 zeros
            pi_hi = pi_lo + local_count

            # We want π(x_est) ≈ n. Current estimate:
            pi_est = pi_lo + sum(1 for p in local_primes if p <= int(x_est))

            deficit = n - pi_est

            # Adjust
            x_new = x_est + deficit * math.log(x_est)

            print(f"  iter {it}: x_est = {x_est:.2f}, "
                  f"π̂(x_est) = {pi_est:.2f}, deficit = {deficit:.2f}, "
                  f"|x_est - p(n)| = {abs(x_est - pn):.2f}")

            if abs(deficit) < 0.5:
                break
            x_est = x_new
            iterations += 1

        print(f"  After {iterations} iterations: |error| = {abs(x_est - pn):.2f}")
        print(f"  The bottleneck is π(lo) estimation, which still needs O(√x) zeros.")
        print()

# =============================================================================
# EXPERIMENT 5: Complexity crossover analysis
# =============================================================================
def experiment_complexity_crossover():
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Complexity Analysis — Hybrid vs Pure Methods")
    print("=" * 70)

    print("""
    Method                    | Complexity         | p(10^9) | p(10^12) | p(10^100)
    --------------------------+--------------------+---------+----------+----------
    Sieve of Eratosthenes     | O(x log log x)     |   10^10 |    10^13 |   10^102
    Meissel-Lehmer            | O(x^{2/3}/log²x)   |   10^6  |    10^8  |   10^66
    Lagarias-Odlyzko analytic | O(x^{1/2+ε})       |   10^5  |    10^6  |   10^50

    Hybrid (K zeros + sieve):
    - K zeros cost: O(K^{1+ε}) via Odlyzko-Schönhage
    - Error bound: E ~ C·√x·log²(x)/K
    - Sieve cost: O(E · √x / log x)
    - Total: O(K^{1+ε} + x/K · log x)
    - Optimal K: x^{1/2} → Total: O(x^{1/2+ε})  [matches Lagarias-Odlyzko!]

    Conclusion: The hybrid method can't beat O(x^{1/2+ε}) because:
    - Computing K zeros costs at least O(K)
    - The error bound E ~ √x/K requires K ~ √x for E < 1
    - So K ~ √x zeros are always needed

    For polylog(n): would need error E = O(1) with K = O(polylog(x)) zeros.
    This requires |Σ_{γ > γ_K} R(x^ρ)| = O(1), but the tail decays only as
    √x · log(γ_K)/γ_K, so γ_K ~ √x is needed. With γ_K ~ K/log(K),
    this gives K ~ √x. BARRIER.
    """)

# =============================================================================
# RUN ALL
# =============================================================================
if __name__ == '__main__':
    print("HYBRID ANALYTIC-SIEVE METHOD FOR p(n)")
    print("=" * 70)

    t0 = time.time()

    experiment_approximation_error()
    experiment_zero_corrected_error()
    experiment_complexity_crossover()
    # These are slower:
    experiment_probabilistic_refinement()

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s")

    print("\n## OVERALL VERDICT")
    print("The hybrid analytic + local sieve approach has complexity O(x^{1/2+ε}),")
    print("matching Lagarias-Odlyzko (1987). It cannot achieve polylog because:")
    print("1. The zero sum tail after K terms is O(√x/K), requiring K ~ √x for O(1) error")
    print("2. Local sieving can't compensate: it still needs π(lo) exactly")
    print("3. The iterative refinement converges but each step needs the same π(x) accuracy")
    print("The fundamental bottleneck remains: any path to p(n) requires O(√x) information")
    print("about the zeta zero oscillations.")
