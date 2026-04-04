"""
Session 6: Prime Zeta Oracle Approach

IDEA: The prime zeta function P(s) = Σ_p p^{-s} encodes ALL primes.
Can we extract individual primes by evaluating P at specific s values?

Key observation: P(s) = Σ_{n=1}^∞ p(n)^{-s}

If we evaluate P(s) for many values of s, we get a system of equations
in the unknowns p(1), p(2), ..., p(N).

But P(s) involves ALL primes, not just the first N. So we'd need to
somehow "truncate" P.

BETTER IDEA: Define P_N(s) = Σ_{n=1}^N p(n)^{-s} (partial prime zeta).
If we know P_N(s) for s = s_1, ..., s_N, can we solve for p(1),...,p(N)?

This is a system of N equations in N unknowns. The unknowns appear as
p(n)^{-s_k}, which is highly nonlinear but has a Vandermonde-like structure.

EVEN BETTER: Work with exponential sums.
Define E(t) = Σ_{n=1}^N e^{-t * p(n)} for t > 0.
If we know E(t) for N values of t, we can in principle recover p(1),...,p(N).
This is the Laplace transform approach / Prony's method.

PRONY'S METHOD: Given a signal f(t) = Σ a_k * e^{-λ_k * t},
recover the frequencies λ_k and amplitudes a_k from samples of f.

For our case: a_k = 1 (all equal), λ_k = p(k).
If we can compute E(t) = Σ e^{-t*p(k)} WITHOUT knowing individual primes,
we can use Prony to extract them.

CAN WE COMPUTE E(t)?
E(t) is related to the partition function of "prime number gas" in statistical mechanics.
Formally: E(t) = Σ_{p prime} e^{-t*p}

Using the explicit formula: Σ_{n≤x} Λ(n) = x - Σ_ρ x^ρ/ρ - ...
The Laplace transform of Λ(n) is -ζ'(s)/ζ(s).
So E(t) is related to evaluations of -ζ'(s)/ζ(s) along certain contours.

This is CIRCULAR again (computing ζ requires knowing primes for the Euler product).
BUT: ζ can be computed from its Dirichlet series for Re(s) > 1, which does NOT
require knowing primes explicitly.
"""

import numpy as np
from mpmath import mp, mpf, zeta, exp, log, pi, sqrt, gamma as mpgamma
from mpmath import nsum, im, re, zetazero
import time

mp.dps = 30

def prime_exp_sum_from_zeta(t, max_x=1000):
    """
    Compute E(t) = Σ_{p prime, p ≤ max_x} e^{-t*p}
    using direct prime enumeration (for testing).
    """
    from sympy import prime, primepi
    N = int(primepi(max_x))
    return sum(float(exp(-t * prime(n))) for n in range(1, N + 1))

def prime_exp_sum_approx(t, method='analytic'):
    """
    Approximate E(t) = Σ_p e^{-tp} using analytic methods.

    By partial summation:
    E(t) = Σ_p e^{-tp} = ∫_2^∞ e^{-tx} dπ(x) = t ∫_2^∞ e^{-tx} π(x) dx + e^{-2t} π(2)

    Using π(x) ≈ li(x):
    E(t) ≈ t ∫_2^∞ e^{-tx} li(x) dx + e^{-2t}

    This integral can be evaluated!
    """
    t = mpf(t)

    if method == 'analytic':
        # E(t) ≈ t ∫_2^∞ e^{-tx}/ln(x) dx + e^{-2t}
        # = t * Ei_1(2t * ln(2)) / ... This is related to exponential integrals

        # Actually: Σ_p p^{-s} = Σ_{k=1}^∞ μ(k)/k * ln(ζ(ks))
        # So E(t) = Σ_p e^{-tp} is the Laplace transform at t

        # For t > 0, we can use the PNT to estimate:
        # E(t) ≈ ∫_2^∞ e^{-tx}/ln(x) dx = Ei(-2t)/(-1) ... not quite

        # Use direct estimate: E(t) ≈ e^{-2t}(1 + e^{-t} + 2*e^{-3t} + ...)
        # = sum of e^{-tp} for small primes + integral for large primes

        result = mpf(0)
        # Small primes explicitly
        small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        for p in small_primes:
            result += exp(-t * p)

        # Large primes: use integral approximation
        # Σ_{p > 47} e^{-tp} ≈ ∫_47^∞ e^{-tx}/ln(x) dx
        from mpmath import quad
        def integrand(x):
            return exp(-t * x) / log(x)

        if t > 0.001:  # Only if integral converges fast enough
            tail = quad(integrand, [mpf(50), mpf(10000)])
            result += tail

        return float(result)

    return 0

def prony_method(samples, t_values, N_components):
    """
    Prony's method: recover frequencies from exponential sum samples.

    Given f(t_k) = Σ_{j=1}^N a_j * exp(-λ_j * t_k),
    recover λ_j and a_j.

    For uniformly spaced t_k = k*dt:
    1. Form the Hankel matrix from samples
    2. Find roots of the characteristic polynomial
    3. Extract frequencies and amplitudes
    """
    M = len(samples)
    if M < 2 * N_components:
        return None

    # Build Hankel matrix
    n = N_components
    H = np.zeros((M - n, n))
    for i in range(M - n):
        for j in range(n):
            H[i, j] = samples[i + j]

    b = samples[n:M]

    # Solve for characteristic polynomial coefficients
    try:
        coeffs, _, _, _ = np.linalg.lstsq(H, b, rcond=None)
    except:
        return None

    # The characteristic polynomial is z^n - c_{n-1}z^{n-1} - ... - c_0
    char_poly = np.zeros(n + 1)
    char_poly[n] = 1
    for i in range(n):
        char_poly[i] = -coeffs[i]

    # Find roots
    roots = np.roots(char_poly[::-1])

    # Convert roots to frequencies: root = exp(-λ * dt)
    dt = t_values[1] - t_values[0]
    frequencies = []
    for r in roots:
        if abs(r) > 1e-10 and r.real > 0:
            freq = -np.log(abs(r)) / dt
            frequencies.append(freq)

    frequencies.sort()
    return frequencies

def test_prony_on_known_primes():
    """Test Prony's method on a known small set of primes."""
    print("="*70)
    print("PRONY'S METHOD: EXTRACTING PRIMES FROM EXPONENTIAL SUMS")
    print("="*70)

    from sympy import prime

    # Generate the "signal" E(t) = Σ_{n=1}^N exp(-t * p(n))
    N_primes = 10  # Try to recover first 10 primes
    true_primes = [int(prime(n)) for n in range(1, N_primes + 1)]
    print(f"\nTrue primes: {true_primes}")

    # Sample E(t) at uniformly spaced points
    n_samples = 3 * N_primes  # Oversample
    dt = 0.05
    t_values = np.array([k * dt for k in range(1, n_samples + 1)])

    # Compute exact samples
    samples = np.array([sum(np.exp(-t * p) for p in true_primes) for t in t_values])

    print(f"\nSamples of E(t) at t = 0.05, 0.10, ..., {t_values[-1]:.2f}:")
    for i in range(min(5, len(samples))):
        print(f"  E({t_values[i]:.2f}) = {samples[i]:.8f}")

    # Apply Prony's method
    recovered = prony_method(samples, t_values, N_primes)
    if recovered:
        print(f"\nRecovered frequencies (should be primes):")
        for i, f in enumerate(recovered):
            nearest_prime = min(true_primes, key=lambda p: abs(p - f))
            print(f"  λ_{i+1} = {f:.4f} (nearest prime: {nearest_prime}, "
                  f"error: {abs(f - nearest_prime):.4f})")

        # Count exact matches
        matched = 0
        for f in recovered:
            for p in true_primes:
                if abs(f - p) < 0.5:
                    matched += 1
                    break
        print(f"\n  Matched: {matched}/{N_primes} primes")
    else:
        print("  Prony's method failed")

    # Try with NOISY samples (as we'd get from analytic approximation)
    print(f"\n--- With noisy samples (1% noise) ---")
    noise = np.random.normal(0, 0.01 * np.abs(samples))
    noisy_samples = samples + noise

    recovered_noisy = prony_method(noisy_samples, t_values, N_primes)
    if recovered_noisy:
        matched_noisy = 0
        for f in recovered_noisy:
            for p in true_primes:
                if abs(f - p) < 0.5:
                    matched_noisy += 1
                    break
        print(f"  Matched with noise: {matched_noisy}/{N_primes} primes")
    else:
        print("  Prony's method failed with noise")

def test_analytic_E():
    """Test the analytic approximation of E(t)."""
    print("\n" + "="*70)
    print("ANALYTIC APPROXIMATION OF E(t)")
    print("="*70)

    for t in [0.01, 0.05, 0.1, 0.5, 1.0]:
        exact = prime_exp_sum_from_zeta(t, max_x=10000)
        approx = prime_exp_sum_approx(t, method='analytic')
        rel_err = abs(exact - approx) / abs(exact) if exact != 0 else float('inf')
        print(f"  t={t:.2f}: exact={exact:.8f}, approx={approx:.8f}, "
              f"rel_err={100*rel_err:.2f}%")

def main():
    print("Session 6: Prime Zeta Oracle / Prony's Method")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    test_prony_on_known_primes()
    test_analytic_E()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
