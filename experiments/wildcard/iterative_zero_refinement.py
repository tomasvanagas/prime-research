#!/usr/bin/env python3
"""
Experiment: Iterative Zero-Sum Refinement for p(n)

NOVEL IDEA: Instead of computing ALL zeta zeros to get π(x) exactly,
use an iterative scheme where:

1. Start: x₀ = Li⁻¹(n) (smooth approximation, ~50% of digits correct)
2. Compute partial zero sum using only T₁ << x zeros
3. Refine: x₁ = Li⁻¹(n + Σ_{|γ|≤T₁} Li(x₀^ρ))
4. The key insight: because the correction depends WEAKLY on x
   (sensitivity ~x^{-1/2}), using x₀ instead of the true p(n)
   introduces only a small error in the correction itself.
5. Iterate with more zeros each round.

QUESTION: How many zeros do we need at each stage?
Can the self-correction property reduce the total zero count?

Also test: Does the partial zero sum converge from below/above?
Is there a way to bound the remaining error after T zeros?
"""

import numpy as np
from sympy import primerange, prime, primepi, li, isprime
from mpmath import mp, mpf, li as mpli, log as mplog, exp as mpexp, pi as mppi, loggamma
import time

# Set precision
mp.dps = 50

# First 100 nontrivial zeros of Riemann zeta (imaginary parts)
# Source: LMFDB / Odlyzko's tables
ZETA_ZEROS = [
    14.134725141734693, 21.022039638771554, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609809, 65.112544048081606,
    67.079810529494173, 69.546401711173979, 72.067157674481907,
    75.704690699083933, 77.144840068874805, 79.337375020249367,
    82.910380854086030, 84.735492980517050, 87.425274613125229,
    88.809111207634465, 92.491899270558484, 94.651344040519838,
    95.870634228245309, 98.831194218193692, 101.31785100573139,
    103.72553804532511, 105.44662305232609, 107.16861118427640,
    111.02953554316967, 111.87465917699263, 114.32022091545271,
    116.22668032085755, 118.79078286597621, 121.37012500242066,
    122.94682929355258, 124.25681855434576, 127.51668387959649,
    129.57870419995605, 131.08768853093265, 133.49773720299758,
    134.75650975337387, 138.11604205453344, 139.73620895212138,
    141.12370740402112, 143.11184580762063,
]

def Li_inverse(n):
    """Compute Li⁻¹(n) using Newton's method on Li(x) = n."""
    # Initial guess from asymptotic: Li⁻¹(n) ≈ n*ln(n)
    x = mpf(n) * mplog(mpf(n))
    for _ in range(100):
        lix = mpli(x)
        if abs(lix - n) < mpf(10)**(-40):
            break
        # Li'(x) = 1/ln(x)
        x = x + (mpf(n) - lix) * mplog(x)
    return x

def partial_zero_sum(x, num_zeros):
    """Compute Σ_{k=1}^{num_zeros} Li(x^ρ_k) + Li(x^{conj(ρ_k)})
    where ρ_k = 1/2 + i*γ_k.

    Using: Li(x^ρ) ≈ x^ρ / (ρ * log(x)) for large x.
    For exact: use the integral representation.
    """
    if num_zeros == 0:
        return mpf(0)

    result = mpf(0)
    logx = mplog(x)

    for k in range(min(num_zeros, len(ZETA_ZEROS))):
        gamma = mpf(ZETA_ZEROS[k])
        rho = mpf(0.5) + 1j * gamma

        # x^ρ = x^{1/2} * x^{iγ} = x^{1/2} * e^{iγ ln x}
        x_half = mpexp(mpf(0.5) * logx)
        phase = gamma * logx
        x_rho = x_half * (mpexp(1j * phase))

        # Li(x^ρ) ≈ Ei(ρ * log x) but we'll use the approximation x^ρ/(ρ log x)
        # which is good for large x
        li_xrho = x_rho / (rho * logx)

        # Add contribution from ρ and conjugate ρ̄
        # Li(x^ρ) + Li(x^ρ̄) = 2*Re(Li(x^ρ))
        result += 2 * li_xrho.real

    return result

def prime_from_zeros(n, max_zeros=50, verbose=True):
    """
    Iterative refinement to compute p(n).

    Step 1: x₀ = Li⁻¹(n)
    Step 2: For increasing T, refine using zero sum
    """
    # Step 1: Smooth approximation
    x0 = Li_inverse(n)
    actual_pn = prime(n)

    if verbose:
        print(f"\n--- Computing p({n}) ---")
        print(f"Actual p({n}) = {actual_pn}")
        print(f"Li⁻¹({n}) = {float(x0):.6f}")
        print(f"Initial error: {float(x0 - actual_pn):.6f}")

    # Step 2: Iterative refinement with increasing zero count
    x_current = x0
    errors = [(0, float(x0 - actual_pn))]

    for num_zeros in [1, 2, 5, 10, 20, 30, 40, 50]:
        if num_zeros > len(ZETA_ZEROS):
            break

        # Compute partial zero sum at current x
        zero_sum = partial_zero_sum(x_current, num_zeros)

        # Refine: x_new = Li⁻¹(n + zero_sum + log(2) - integral_term)
        # Simplified: π(x) ≈ Li(x) - zero_sum - log(2)
        # So: Li(p(n)) ≈ n + zero_sum + log(2)
        correction = zero_sum + mplog(2)
        x_new = Li_inverse(float(mpf(n) + correction))

        error = float(x_new - actual_pn)
        errors.append((num_zeros, error))

        if verbose:
            print(f"  {num_zeros:3d} zeros: x = {float(x_new):.6f}, error = {error:.6f}")

        x_current = x_new

    # Step 3: Self-correction iteration (use x_new to recompute zero sum)
    if verbose:
        print(f"\n  Self-correction iterations (50 zeros):")

    for iteration in range(5):
        zero_sum = partial_zero_sum(x_current, min(50, len(ZETA_ZEROS)))
        correction = zero_sum + mplog(2)
        x_new = Li_inverse(float(mpf(n) + correction))
        error = float(x_new - actual_pn)

        if verbose:
            print(f"    iter {iteration}: error = {error:.6f}")

        x_current = x_new

    return float(x_current), errors

def analyze_zero_convergence():
    """
    KEY QUESTION: How does the error decrease as we add more zeros?
    If it decreases faster than 1/T, we might beat Lagarias-Odlyzko.
    """
    print("=" * 70)
    print("ZERO SUM CONVERGENCE ANALYSIS")
    print("=" * 70)

    test_cases = [100, 500, 1000, 5000, 10000]

    for n in test_cases:
        actual_pn = prime(n)
        x0 = Li_inverse(n)

        print(f"\nn = {n}, p(n) = {actual_pn}")
        print(f"  Li⁻¹(n) = {float(x0):.2f}, error = {float(x0 - actual_pn):.2f}")

        # Compute zero sum at TRUE p(n) with increasing zeros
        errors = []
        for nz in range(0, min(51, len(ZETA_ZEROS) + 1)):
            zs = partial_zero_sum(mpf(actual_pn), nz)
            # π(x) = Li(x) - zero_sum - log(2) + ...
            # So Li(x) - π(x) should approach zero_sum + log(2)
            estimated_pi = float(mpli(mpf(actual_pn)) - zs - mplog(2))
            error = estimated_pi - n  # Should approach 0
            errors.append((nz, error))

        print(f"  Errors with increasing zeros (at true p(n)):")
        for nz, err in errors[::5]:
            print(f"    {nz:3d} zeros: π_approx - n = {err:.4f}")

        # What's the minimum number of zeros for |error| < 0.5?
        for nz, err in errors:
            if abs(err) < 0.5:
                print(f"  → Exact with {nz} zeros (error {err:.4f} < 0.5)")
                break
        else:
            print(f"  → NOT exact with {len(ZETA_ZEROS)} zeros (last error: {errors[-1][1]:.4f})")

def test_sensitivity():
    """
    Test: How sensitive is the zero sum to the evaluation point?
    If ∂(zero_sum)/∂x is small, self-correction converges fast.
    """
    print("\n" + "=" * 70)
    print("SENSITIVITY ANALYSIS")
    print("=" * 70)

    for n in [100, 1000, 10000]:
        actual_pn = prime(n)
        x0 = Li_inverse(n)
        delta = float(x0 - actual_pn)

        # Compute zero sum at p(n) and at Li⁻¹(n)
        zs_true = partial_zero_sum(mpf(actual_pn), 50)
        zs_approx = partial_zero_sum(x0, 50)
        diff = float(zs_approx - zs_true)

        print(f"\nn={n}, p(n)={actual_pn}, Li⁻¹(n)={float(x0):.2f}")
        print(f"  δ = Li⁻¹(n) - p(n) = {delta:.4f}")
        print(f"  Zero sum at p(n):    {float(zs_true):.6f}")
        print(f"  Zero sum at Li⁻¹(n): {float(zs_approx):.6f}")
        print(f"  Difference: {diff:.6f}")
        print(f"  Sensitivity: Δ(zero_sum)/δ = {diff/delta:.6f} (should be << 1)")
        print(f"  Relative: |diff/zero_sum| = {abs(diff/float(zs_true)) if float(zs_true) != 0 else 'inf':.6f}")

def novel_contour_approach():
    """
    NOVEL APPROACH: Instead of summing individual zeros,
    evaluate the contour integral directly.

    π(x) = (1/2πi) ∮ log(ζ(s)) x^s / s ds

    Deform the contour to the right of Re(s) = 1 where ζ(s) is
    computed via Euler product. The integral then involves only
    primes, not zeros.

    But we need to separate the contribution from primes ≤ x
    (which is what we want) from primes > x.
    """
    print("\n" + "=" * 70)
    print("CONTOUR INTEGRAL APPROACH")
    print("=" * 70)

    # Evaluate: (1/2πi) ∫_{c-iT}^{c+iT} [-ζ'/ζ(s)] x^s/s ds
    # For c > 1, ζ'/ζ(s) = -Σ_n Λ(n)/n^s converges.
    # The integral picks up Ψ(x) = Σ_{n≤x} Λ(n)

    # Test: compute Ψ(x) via numerical contour integration
    from sympy import primepi

    for x in [50, 100, 500]:
        # True Ψ(x)
        true_psi = sum(np.log(p) * int(np.log(x) / np.log(p))
                       for p in primerange(2, x+1))

        # Numerical contour integral
        c = 2.0  # Contour at Re(s) = c
        T = 100.0  # Height of contour
        N_quad = 1000  # Quadrature points

        ts = np.linspace(-T, T, N_quad)
        dt = ts[1] - ts[0]

        # Evaluate -ζ'/ζ at s = c + it using Euler product (truncated)
        integral = 0.0
        primes_for_euler = list(primerange(2, 1000))

        for t in ts:
            s = complex(c, t)
            # -ζ'/ζ(s) ≈ Σ_p log(p)/(p^s - 1) (from Euler product)
            zeta_log_deriv = sum(np.log(p) / (p**s - 1) for p in primes_for_euler)

            # x^s / s
            integrand = zeta_log_deriv * x**s / s

            integral += integrand * dt

        psi_approx = (integral / (2 * np.pi * 1j)).real
        error = psi_approx - true_psi

        print(f"\nx={x}: Ψ(x) = {true_psi:.4f}")
        print(f"  Contour integral (c={c}, T={T}, N={N_quad}): {psi_approx:.4f}")
        print(f"  Error: {error:.4f}")
        print(f"  Relative error: {abs(error/true_psi)*100:.2f}%")

def test_iterative_full():
    """Full test of iterative refinement on various n."""
    print("\n" + "=" * 70)
    print("FULL ITERATIVE REFINEMENT TEST")
    print("=" * 70)

    results = []
    for n in [10, 50, 100, 500, 1000, 5000, 10000]:
        x_final, errors = prime_from_zeros(n, verbose=False)
        actual = prime(n)
        final_error = x_final - actual

        # Find minimum zeros for |error| < 1
        min_zeros = None
        for nz, err in errors:
            if abs(err) < 1.0:
                min_zeros = nz
                break

        results.append((n, actual, float(Li_inverse(n)), x_final, final_error, min_zeros))
        print(f"n={n:6d}: p(n)={actual:8d}, Li⁻¹={float(Li_inverse(n)):10.2f}, "
              f"refined={x_final:10.2f}, error={final_error:+8.3f}, "
              f"zeros_needed={min_zeros}")

    # Analyze: does zeros_needed grow with n? How fast?
    print("\nScaling of zeros needed:")
    for n, _, _, _, _, mz in results:
        if mz is not None:
            print(f"  n={n}: {mz} zeros (√p(n) ≈ {int(prime(n)**0.5)})")


if __name__ == "__main__":
    print("ITERATIVE ZERO-SUM REFINEMENT FOR p(n)")
    print("Testing if partial zero sums + self-correction = exact p(n)\n")

    # Test 1: Basic iterative refinement
    for n in [100, 1000]:
        prime_from_zeros(n, verbose=True)

    # Test 2: Convergence analysis
    analyze_zero_convergence()

    # Test 3: Sensitivity (key to self-correction)
    test_sensitivity()

    # Test 4: Contour integral approach
    novel_contour_approach()

    # Test 5: Full test
    test_iterative_full()

    print("\n" + "=" * 70)
    print("KEY QUESTIONS ANSWERED")
    print("=" * 70)
    print("""
    1. How many zeros for exact p(n)?
       → If it's O(polylog(n)), we have a breakthrough!
       → If it's O(√n) or O(n^ε), it's the known barrier.

    2. Does self-correction converge?
       → If sensitivity < 1, iteration converges.
       → Rate of convergence determines practical speed.

    3. Can the contour integral bypass individual zeros?
       → If quadrature points << number of zeros, this wins.
    """)
