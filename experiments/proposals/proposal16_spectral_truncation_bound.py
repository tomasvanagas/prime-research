#!/usr/bin/env python3
"""
PROPOSAL 16: Spectral Truncation with Adaptive Zero Selection
═══════════════════════════════════════════════════════════════

IDEA: Instead of summing over ALL zeta zeros in the explicit formula,
select a small O(polylog(n)) subset of zeros that contribute the most
to delta(n) = p(n) - R^{-1}(n). Use the GUE statistics of zero spacing
to predict which zeros matter most for a given n, then compute only those
contributions.

MATHEMATICAL BASIS:
  pi(x) = R(x) - sum_{rho} R(x^rho)

  The contribution of zero rho = 1/2 + i*gamma to R(x^rho) decays like
  x^{1/2} / |rho| * cos(gamma * log(x) + phase).

  For a specific x = p(n), the oscillatory sum has cancellations.
  Key insight: zeros whose gamma * log(x) is near a multiple of 2*pi
  contribute coherently. We can identify these "resonant zeros" in
  O(polylog) time using continued fraction approximation.

CONJECTURE: For most n, O(log^k(n)) carefully selected zeros suffice
to determine delta(n) exactly (or to within < 0.5, allowing rounding).

TEST: For small n, compare the explicit formula truncated at various
subsets of zeros against the true delta(n).
"""

import math
import sys
from functools import lru_cache

# Load zeta zeros
def load_zeros(filepath):
    zeros = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        pass
    return zeros

# Riemann R function via Gram series: R(x) = sum_{k=0}^{inf} (ln x)^k / (k! * k * zeta(k+1))
# with the convention that the k=0 term is 1. Actually:
# R(x) = 1 + sum_{k=1}^{inf} (ln x)^k / (k * k! * zeta(k+1))
def riemann_R(x):
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    total = 1.0
    log_power = 1.0  # (lnx)^k / k!
    for k in range(1, 200):
        log_power *= lnx / k
        zk = _zeta_approx(k + 1)
        term = log_power / (k * zk)
        total += term
        if abs(term) < 1e-15:
            break
    return total

def _zeta_approx(s):
    """Approximate zeta(s) for s >= 2."""
    if s == 2:
        return math.pi**2 / 6
    if s == 3:
        return 1.2020569031595942
    if s == 4:
        return math.pi**4 / 90
    # For large s, zeta(s) ~ 1
    total = 0.0
    for k in range(1, 100):
        total += 1.0 / k**s
        if 1.0 / k**s < 1e-15:
            break
    return total

# Inverse Riemann R: given n, find x such that R(x) = n
def riemann_R_inv(n):
    if n <= 1:
        return 2.0
    # Initial guess: n * log(n)
    x = n * math.log(n) + n * math.log(math.log(n + 2))
    for _ in range(200):
        rx = riemann_R(x)
        if abs(rx - n) < 0.0001:
            break
        # Newton step: R'(x) ~ 1/log(x)
        deriv = 1.0 / math.log(x) if x > 1 else 1.0
        x += (n - rx) / deriv
        x = max(x, 2.0)
    return x

# True nth prime (sieve for small n)
def true_nth_prime(n):
    if n < 1:
        return 2
    limit = max(100, int(n * (math.log(n) + math.log(math.log(n + 2)) + 3)))
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    primes = [i for i in range(2, limit + 1) if sieve[i]]
    if n <= len(primes):
        return primes[n - 1]
    return None

# Contribution of a single zero to R(x^rho)
def zero_contribution(gamma, x):
    """
    R(x^rho) where rho = 1/2 + i*gamma.
    Leading term: x^{1/2} * cos(gamma * ln(x)) / |rho| approximately.
    More precisely, use the series expansion.
    """
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    sqrt_x = math.sqrt(x)
    # Leading oscillatory contribution
    phase = gamma * lnx
    # The contribution magnitude decays as 1/|rho|
    magnitude = sqrt_x / math.sqrt(0.25 + gamma**2)
    return 2.0 * magnitude * math.cos(phase)  # Factor 2 for conjugate pair

# Resonance score: how coherently does this zero contribute at x?
def resonance_score(gamma, x):
    """
    Zeros are most impactful when gamma * log(x) is near k*pi for integer k.
    The contribution is proportional to |cos(gamma * log(x))| / |rho|.
    """
    phase = gamma * math.log(x)
    return abs(math.cos(phase)) / math.sqrt(0.25 + gamma**2)

# Select top-K resonant zeros
def select_resonant_zeros(zeros, x, K):
    scored = [(resonance_score(g, x), g) for g in zeros]
    scored.sort(reverse=True)
    return [g for _, g in scored[:K]]

# Compute delta(n) using selected zeros
def compute_delta_from_zeros(n, zeros_subset, x_estimate):
    """Sum zero contributions to estimate correction."""
    correction = sum(zero_contribution(g, x_estimate) for g in zeros_subset)
    return correction

def run_test():
    # Load available zeros
    zeros = load_zeros("data/zeta_zeros_1000.txt")
    if not zeros:
        zeros = load_zeros("../../data/zeta_zeros_1000.txt")
    if not zeros:
        # Use first 50 known zeros
        zeros = [
            14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
            37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
            52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
            67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
            79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
            92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
            103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
            114.320220, 116.226680, 118.790783, 121.370125, 122.946829,
            124.256819, 127.516684, 129.578704, 131.087688, 133.497737,
            134.756510, 138.116042, 139.736209, 141.123707, 143.111846,
        ]

    print("PROPOSAL 16: Spectral Truncation with Adaptive Zero Selection")
    print("=" * 65)
    print(f"Available zeros: {len(zeros)}")
    print()

    # Test for various n
    test_ns = [10, 50, 100, 500, 1000, 2000, 5000, 10000]
    K_values = [5, 10, 20, 50]

    results = {}

    for n in test_ns:
        true_p = true_nth_prime(n)
        if true_p is None:
            continue

        x_est = riemann_R_inv(n)
        delta_true = true_p - x_est

        results[n] = {
            'true_p': true_p,
            'R_inv': x_est,
            'delta_true': delta_true,
            'K_results': {}
        }

        for K in K_values:
            if K > len(zeros):
                continue
            selected = select_resonant_zeros(zeros, x_est, K)
            delta_est = compute_delta_from_zeros(n, selected, x_est)

            # Also try ALL first K zeros (non-adaptive)
            delta_sequential = compute_delta_from_zeros(n, zeros[:K], x_est)

            error_adaptive = abs(delta_true - delta_est)
            error_sequential = abs(delta_true - delta_sequential)

            results[n]['K_results'][K] = {
                'delta_adaptive': delta_est,
                'delta_sequential': delta_sequential,
                'error_adaptive': error_adaptive,
                'error_sequential': error_sequential,
                'correct_adaptive': abs(round(x_est + delta_est) - true_p) == 0,
                'correct_sequential': abs(round(x_est + delta_sequential) - true_p) == 0,
            }

    # Print results table
    print(f"{'n':>6} | {'p(n)':>8} | {'R_inv':>10} | {'delta':>8} | K=5 err | K=10 err | K=20 err | K=50 err")
    print("-" * 90)

    adaptive_wins = 0
    total_tests = 0
    exact_with_few = 0

    for n in test_ns:
        if n not in results:
            continue
        r = results[n]
        row = f"{n:>6} | {r['true_p']:>8} | {r['R_inv']:>10.2f} | {r['delta_true']:>8.2f} |"
        for K in K_values:
            if K in r['K_results']:
                kr = r['K_results'][K]
                row += f" {kr['error_adaptive']:>7.2f} |"
                if kr['error_adaptive'] < kr['error_sequential']:
                    adaptive_wins += 1
                if kr['correct_adaptive'] and K <= 20:
                    exact_with_few += 1
                total_tests += 1
            else:
                row += "     N/A |"
        print(row)

    print()
    print(f"Adaptive selection beats sequential: {adaptive_wins}/{total_tests}")
    print(f"Exact with K<=20 zeros (adaptive): {exact_with_few} cases")

    # Detailed analysis: how many zeros needed for exactness?
    print()
    print("ZEROS NEEDED FOR EXACT p(n):")
    print(f"{'n':>6} | {'zeros_needed':>12} | {'log(n)':>8} | {'log^2(n)':>8}")
    print("-" * 50)

    for n in test_ns:
        if n not in results:
            continue
        r = results[n]
        x_est = r['R_inv']

        # Find minimum K for exact answer (adaptive)
        min_K = None
        for K in range(1, min(len(zeros) + 1, 200)):
            selected = select_resonant_zeros(zeros, x_est, K)
            delta_est = compute_delta_from_zeros(n, selected, x_est)
            if round(x_est + delta_est) == r['true_p']:
                min_K = K
                break

        logn = math.log(n) if n > 1 else 1
        log2n = logn ** 2

        k_str = str(min_K) if min_K else ">all"
        print(f"{n:>6} | {k_str:>12} | {logn:>8.2f} | {log2n:>8.2f}")

    return results

if __name__ == "__main__":
    results = run_test()
