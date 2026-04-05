#!/usr/bin/env python3
"""
Proposal 21: Zero Clustering Truncation for Explicit Formula

IDEA: The Riemann explicit formula for pi(x) requires summing over ALL nontrivial
zeros rho = 1/2 + i*gamma of zeta(s). Normally this needs O(sqrt(x)) zeros.

Key insight: The zeros exhibit GUE (Gaussian Unitary Ensemble) statistics, meaning
their local spacings follow predictable distributions. We can exploit this to:

1. Group zeros into "clusters" where nearby zeros have correlated contributions
2. Represent each cluster by a single representative term + correction
3. Use the pair correlation function R_2(u) = 1 - (sin(pi*u)/(pi*u))^2 to
   predict the aggregate contribution of unseen zeros from seen ones

The hope: If cluster contributions decay fast enough, we might need only
O(polylog(x)) representative clusters instead of O(sqrt(x)) individual zeros.

CONJECTURE: The oscillatory sum S(x) = sum_gamma x^(i*gamma)/rho converges
to within O(1) accuracy using O(log(x)^C) carefully chosen zeros, for some C.

TIME COMPLEXITY: O(log(x)^C * M(log x)) where M is multiplication cost.
"""

import math
import os
from functools import lru_cache

# Load zeta zeros
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')

def load_zeros(n=200):
    path = os.path.join(DATA_DIR, f'zeta_zeros_{n}.txt')
    zeros = []
    with open(path) as f:
        for line in f:
            val = line.strip()
            if val:
                zeros.append(float(val))
    return zeros

ZEROS = load_zeros(200)

def R_inverse(n):
    """Approximate inverse of Riemann R function using Newton's method."""
    # Initial estimate from PNT
    x = n * math.log(n) + n * math.log(math.log(n)) - n
    if x < 2:
        x = 2.0
    for _ in range(50):
        rx = R_function(x)
        rx_prime = R_function_derivative(x)
        if abs(rx_prime) < 1e-15:
            break
        x_new = x - (rx - n) / rx_prime
        if abs(x_new - x) < 1e-10:
            break
        x = x_new
    return x

def R_function(x):
    """Riemann R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})"""
    if x <= 1:
        return 0.0
    total = 0.0
    # Mobius function values for k=1..20
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    for k in range(1, 21):
        if mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            break
        total += mu[k] / k * li(xk)
    return total

def R_function_derivative(x):
    """Derivative of R(x) with respect to x."""
    if x <= 1:
        return 0.0
    total = 0.0
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0]
    for k in range(1, 21):
        if mu[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001:
            break
        # d/dx [li(x^{1/k})] = 1/(k * log(x^{1/k})) * x^{1/k - 1}
        #                     = 1/(log(x)) * x^{-1}  (simplified)
        # Actually: d/dx li(x^{1/k}) = x^{1/k - 1} / (k * ln(x^{1/k}))
        logxk = math.log(xk)
        if abs(logxk) < 1e-15:
            continue
        total += mu[k] / k * (xk / x) / (k * logxk)
    return total

def li(x):
    """Logarithmic integral li(x) = integral from 0 to x of dt/ln(t)."""
    if x <= 0:
        return 0.0
    if x <= 1:
        return 0.0
    if abs(x - 1.0) < 0.01:
        return -100.0  # singularity
    # Numerical approximation using Ramanujan's series
    lnx = math.log(x)
    total = 0.0
    term = 1.0
    for k in range(1, 100):
        term *= lnx / k
        total += term / k
        if abs(term / k) < 1e-15:
            break
    return 0.5772156649015329 + math.log(abs(lnx)) + total

def explicit_formula_correction(x, num_zeros):
    """
    Compute oscillatory correction from Riemann explicit formula.
    pi(x) ≈ R(x) - sum_{rho} R(x^rho)

    For each zero gamma, the contribution is approximately:
    -2 * Re(li(x^{1/2 + i*gamma})) / 1  (simplified)

    Using the approximation: li(x^rho) ≈ x^rho / (rho * log(x))
    """
    if x <= 2:
        return 0.0
    lnx = math.log(x)
    correction = 0.0
    for j in range(min(num_zeros, len(ZEROS))):
        gamma = ZEROS[j]
        # x^{i*gamma} = e^{i*gamma*ln(x)} = cos(gamma*ln(x)) + i*sin(gamma*ln(x))
        # x^{1/2 + i*gamma} = sqrt(x) * x^{i*gamma}
        # Re[li(x^rho)] ≈ Re[x^rho / (rho * ln(x))]
        # rho = 1/2 + i*gamma, |rho|^2 = 1/4 + gamma^2
        # Re[x^rho / rho] = sqrt(x) * (cos(gamma*lnx)/2 + gamma*sin(gamma*lnx)) / (1/4 + gamma^2)
        phase = gamma * lnx
        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)
        sqrtx = math.sqrt(x)
        rho_norm_sq = 0.25 + gamma * gamma
        # Contribution from rho and conj(rho) combined: -2 * Re[li(x^rho)]
        real_part = sqrtx * (0.5 * cos_phase + gamma * sin_phase) / rho_norm_sq
        correction -= 2.0 * real_part / lnx
    return correction

def cluster_zeros(zeros, cluster_radius=2.0):
    """
    Group zeros into clusters based on proximity.
    Returns list of (center, weight, spread) tuples.
    """
    clusters = []
    i = 0
    while i < len(zeros):
        cluster_start = i
        cluster_sum = zeros[i]
        cluster_count = 1
        while i + 1 < len(zeros) and zeros[i + 1] - zeros[cluster_start] < cluster_radius:
            i += 1
            cluster_sum += zeros[i]
            cluster_count += 1
        center = cluster_sum / cluster_count
        spread = zeros[i] - zeros[cluster_start] if cluster_count > 1 else 0
        clusters.append((center, cluster_count, spread))
        i += 1
    return clusters

def clustered_correction(x, clusters):
    """
    Compute correction using clustered zeros.
    Each cluster contributes: weight * individual_contribution(center)
    with a spread correction factor.
    """
    if x <= 2:
        return 0.0
    lnx = math.log(x)
    sqrtx = math.sqrt(x)
    correction = 0.0
    for center, weight, spread in clusters:
        phase = center * lnx
        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)
        rho_norm_sq = 0.25 + center * center
        real_part = sqrtx * (0.5 * cos_phase + center * sin_phase) / rho_norm_sq
        # Spread correction: for a cluster of width w, the phases spread by w*lnx
        # If w*lnx is small, contributions add coherently (weight factor)
        # If w*lnx is large, they partially cancel (sinc-like factor)
        phase_spread = spread * lnx / 2.0
        if phase_spread > 0.01:
            coherence = math.sin(phase_spread) / phase_spread
        else:
            coherence = 1.0
        correction -= 2.0 * weight * coherence * real_part / lnx
    return correction

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(limit + 1) if sieve[i]]

def nth_prime_explicit(n, num_zeros=200):
    """Compute p(n) using R^{-1}(n) + explicit formula correction."""
    x_est = R_inverse(n)
    # Use explicit formula to refine
    correction = explicit_formula_correction(x_est, num_zeros)
    pi_est = R_function(x_est) + correction
    # Binary search to find exact p(n)
    # The correction tells us approximately how far pi(x_est) is from n
    return x_est, pi_est

def test_accuracy(max_n=1000):
    """Test accuracy of various truncation levels."""
    primes = sieve_primes(max_n * 15)  # generous upper bound

    print("=" * 80)
    print("PROPOSAL 21: Zero Clustering Truncation")
    print("=" * 80)

    # Test 1: How many zeros needed for given accuracy?
    print("\n--- Test 1: Accuracy vs number of zeros ---")
    print(f"{'n':>6} {'p(n)':>8} {'R_inv':>10} {'delta':>8} ", end="")
    for nz in [5, 10, 20, 50, 100, 200]:
        print(f"{'corr_'+str(nz):>10}", end="")
    print()

    test_ns = [10, 50, 100, 500, 1000, 2000, 5000]
    for n in test_ns:
        if n > len(primes):
            break
        pn = primes[n - 1]
        x_est = R_inverse(n)
        delta = pn - x_est
        print(f"{n:>6} {pn:>8} {x_est:>10.2f} {delta:>8.2f} ", end="")
        for nz in [5, 10, 20, 50, 100, 200]:
            corr = explicit_formula_correction(x_est, nz)
            # pi(x_est) ≈ R(x_est) + correction ≈ n + correction
            # So the correction to x is approximately correction * ln(pn)
            print(f"{corr:>10.2f}", end="")
        print()

    # Test 2: Clustered vs individual zeros
    print("\n--- Test 2: Clustered zeros vs individual zeros ---")
    for radius in [1.0, 2.0, 5.0, 10.0]:
        clusters = cluster_zeros(ZEROS, cluster_radius=radius)
        print(f"\nCluster radius={radius}: {len(clusters)} clusters from {len(ZEROS)} zeros")

        errors_individual = []
        errors_clustered = []
        for n in [100, 500, 1000, 5000]:
            if n > len(primes):
                break
            pn = primes[n - 1]
            x_est = R_inverse(n)

            corr_ind = explicit_formula_correction(x_est, len(ZEROS))
            corr_cls = clustered_correction(x_est, clusters)

            pi_ind = R_function(x_est) + corr_ind
            pi_cls = R_function(x_est) + corr_cls

            err_ind = abs(pi_ind - n)
            err_cls = abs(pi_cls - n)
            errors_individual.append(err_ind)
            errors_clustered.append(err_cls)
            print(f"  n={n:>5}: individual err={err_ind:.3f}, clustered err={err_cls:.3f}")

    # Test 3: Can we predict the contribution of zeros > T from zeros < T?
    print("\n--- Test 3: Tail prediction via GUE statistics ---")
    # The mean density of zeros at height T is (1/2pi)*log(T/2pi)
    # The pair correlation function gives us information about correlations
    print("Testing whether low zeros can predict high-zero contributions...")

    for n in [100, 500, 1000]:
        if n > len(primes):
            break
        pn = primes[n - 1]
        x_est = R_inverse(n)

        # Full correction with all zeros
        full_corr = explicit_formula_correction(x_est, 200)

        # Partial corrections
        for cutoff in [10, 20, 50]:
            partial = explicit_formula_correction(x_est, cutoff)
            tail = full_corr - partial
            # Try to predict tail from partial info
            # Simple model: tail contribution decays as sum 1/gamma for gamma > cutoff_height
            cutoff_height = ZEROS[cutoff - 1] if cutoff <= len(ZEROS) else ZEROS[-1]
            remaining_density = math.log(cutoff_height / (2 * math.pi)) / (2 * math.pi)
            # This is a very rough model
            predicted_tail = 0.0  # placeholder for sophisticated prediction

            print(f"  n={n}, {cutoff} zeros: partial={partial:.3f}, "
                  f"tail={tail:.3f}, ratio=tail/partial={tail/partial:.3f}" if abs(partial) > 0.001
                  else f"  n={n}, {cutoff} zeros: partial={partial:.3f}, tail={tail:.3f}")

    # Test 4: Error statistics across many n values
    print("\n--- Test 4: Error distribution for R^{-1}(n) vs explicit formula ---")
    errors_rinv = []
    errors_10z = []
    errors_50z = []
    errors_200z = []

    for n in range(10, min(5001, len(primes) + 1)):
        pn = primes[n - 1]
        x_est = R_inverse(n)
        errors_rinv.append(abs(pn - x_est))

        pi_10 = R_function(x_est) + explicit_formula_correction(x_est, 10)
        pi_50 = R_function(x_est) + explicit_formula_correction(x_est, 50)
        pi_200 = R_function(x_est) + explicit_formula_correction(x_est, 200)

        errors_10z.append(abs(pi_10 - n))
        errors_50z.append(abs(pi_50 - n))
        errors_200z.append(abs(pi_200 - n))

    import statistics
    for label, errs in [("R^{-1} only", errors_rinv),
                         ("10 zeros", errors_10z),
                         ("50 zeros", errors_50z),
                         ("200 zeros", errors_200z)]:
        print(f"  {label:>15}: mean={statistics.mean(errs):.3f}, "
              f"median={statistics.median(errs):.3f}, "
              f"max={max(errs):.3f}, "
              f"std={statistics.stdev(errs):.3f}")

    print("\n--- VERDICT ---")
    print("If clustered corrections match individual to within O(1),")
    print("and if cluster count grows as O(log(x)^C), this could work.")
    print("Otherwise the approach reduces constants but not asymptotics.")

if __name__ == '__main__':
    test_accuracy()
