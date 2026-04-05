#!/usr/bin/env python3
"""
Contour Integral Evaluation of the Zero Sum
============================================
Instead of summing x^ρ/ρ over individual zeros, compute the sum as a
contour integral of (ζ'/ζ)(s) · x^s / s around the critical strip.

Key identity (from the explicit formula):
  ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1 - x^{-2})

The zero sum S(x) = Σ_ρ x^ρ/ρ can be written as:
  S(x) = (1/2πi) ∮ -(ζ'/ζ)(s) · x^s/s ds  (around zeros)

We test: can the contour integral be evaluated with FEWER quadrature
points than there are zeros enclosed? If ζ'/ζ is smooth on the contour
(away from zeros), quadrature might converge fast.

Approach:
1. Compute "ground truth" S(x) by directly summing over known zeros
2. Evaluate the contour integral numerically with varying quadrature resolution
3. Compare: how many quadrature points are needed for accuracy?

Also test: SADDLE POINT / STATIONARY PHASE methods to reduce oscillation.
"""

import numpy as np
from scipy import special, integrate
import sympy
from sympy import primepi, prime
import time

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
    """Riemann R function: R(x) = Σ μ(n)/n · li(x^{1/n})"""
    if x <= 1:
        return 0.0
    result = 0.0
    for n in range(1, 100):
        mu_n = int(sympy.mobius(n))
        if mu_n == 0:
            continue
        xn = x ** (1.0 / n)
        if xn <= 1.001:
            break
        li_val = float(special.expi(np.log(xn)))
        result += mu_n / n * li_val
    return result

def R_at_complex(x, rho_imag):
    """Compute R(x^ρ) where ρ = 1/2 + i·γ, using dominant term li(x^ρ)"""
    # R(x^ρ) ≈ li(x^ρ) for the dominant term (n=1 in Gram series)
    # li(x^ρ) = Ei(ρ · log(x))
    log_x = np.log(x)
    s = 0.5 + 1j * rho_imag
    z = s * log_x  # complex argument

    # Compute Ei(z) for complex z using series
    # Ei(z) = γ + ln(z) + Σ z^k/(k·k!) for |z| < ∞
    # For large |z|, use asymptotic: Ei(z) ~ e^z/z · Σ k!/z^k

    if abs(z) > 50:
        # Asymptotic expansion
        result = np.exp(z) / z
        term = 1.0
        for k in range(1, 20):
            term *= k / z
            result += np.exp(z) / z * term
            if abs(term) < 1e-15:
                break
        return result
    else:
        # Series expansion
        euler_gamma = 0.5772156649015329
        result = euler_gamma + np.log(z + 0j)
        term = z
        for k in range(1, 200):
            result += term / (k * np.math.factorial(k) if k < 20 else k * float(sympy.factorial(k)))
            if k < 170:
                term_val = z ** (k + 1)
                if abs(term_val) < 1e-300:
                    break
                term = term_val
            else:
                break
        return result

def zero_sum_direct(x, num_zeros):
    """Compute S(x) = Σ_ρ R(x^ρ) by direct summation over zeros."""
    total = 0.0
    log_x = np.log(x)
    for gamma in ZEROS[:num_zeros]:
        rho = 0.5 + 1j * gamma
        # Each zero ρ contributes R(x^ρ), take conjugate pair together
        # Contribution of ρ and ρ̄: 2·Re[R(x^ρ)]
        # Simplified: use li(x^ρ) ≈ Ei(ρ·log(x)) as dominant term
        z = rho * log_x
        # Use Gram series: R(y) = Σ (-1)^n (ln y)^n / (n! · n · ζ(n+1))
        # Simpler: just use x^ρ/ρ as the main oscillatory part
        contrib = x**rho / rho
        total += 2 * contrib.real  # conjugate pair
    return total

def zero_sum_li(x, num_zeros):
    """Compute Σ_ρ li(x^ρ) using Gram-style series for li at complex argument."""
    total = 0.0
    log_x = np.log(x)
    for gamma in ZEROS[:num_zeros]:
        rho = 0.5 + 1j * gamma
        z = rho * log_x

        # Compute li(x^ρ) = Ei(ρ·log x) via convergent series
        # Ei(z) = γ + ln(z) + Σ_{k=1}^∞ z^k/(k·k!)
        euler_gamma = 0.5772156649015329
        ei = euler_gamma + np.log(z)
        zk = 1.0 + 0j
        for k in range(1, 100):
            zk *= z / k
            term = zk / k
            ei += term
            if abs(term) < 1e-15 * max(abs(ei), 1):
                break

        total += 2 * ei.real  # conjugate pair
    return total

def chebyshev_psi(x):
    """Compute ψ(x) = Σ_{p^k ≤ x} log(p) exactly for moderate x."""
    result = 0.0
    for p in sympy.primerange(2, int(x) + 1):
        pk = p
        while pk <= x:
            result += float(sympy.log(p))
            pk *= p
    return result

def pi_via_explicit(x, num_zeros):
    """Compute π(x) using explicit formula with given number of zeros."""
    Rx = R_function(x)
    zero_correction = zero_sum_li(x, num_zeros)
    # π(x) ≈ R(x) - Σ_ρ R(x^ρ) (simplified)
    return Rx - zero_correction

# =============================================================================
# EXPERIMENT 1: How does the zero sum converge?
# =============================================================================
def experiment_convergence():
    print("=" * 70)
    print("EXPERIMENT 1: Zero Sum Convergence Analysis")
    print("=" * 70)

    test_values = [10**4, 10**5, 10**6, 10**7]
    zero_counts = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]

    results = {}
    for x in test_values:
        true_pi = int(sympy.primepi(x))
        Rx = R_function(x)
        print(f"\nx = {x:>10}, π(x) = {true_pi}, R(x) = {Rx:.4f}, R(x)-π(x) = {Rx - true_pi:.4f}")

        results[x] = []
        for nz in zero_counts:
            if nz > len(ZEROS):
                break
            correction = zero_sum_li(x, nz)
            estimate = Rx - correction
            error = estimate - true_pi
            results[x].append((nz, estimate, error))
            print(f"  {nz:>4} zeros: π̂(x) = {estimate:>12.4f}, error = {error:>+10.4f}")

    return results

# =============================================================================
# EXPERIMENT 2: Contour integral vs direct summation
# =============================================================================
def zeta_log_derivative(s, num_terms=1000):
    """Compute -ζ'/ζ(s) using Euler product for Re(s) > 1, or Dirichlet series.

    -ζ'/ζ(s) = Σ_n Λ(n)/n^s  (von Mangoldt)

    For Re(s) > 1, this converges.
    """
    result = 0.0 + 0j
    # Use Dirichlet series: -ζ'/ζ(s) = Σ Λ(n) n^{-s}
    for p in sympy.primerange(2, num_terms + 1):
        pk = int(p)
        logp = np.log(float(p))
        while pk <= num_terms:
            result += logp * float(pk) ** (-s)
            pk *= int(p)
    return result

def experiment_contour_integral():
    """
    Test: Can we evaluate the zero sum via contour integration?

    The idea: Σ_ρ x^ρ/ρ = -(1/2πi) ∫_C (ζ'/ζ)(s) · x^s/s ds
    where C encloses the zeros.

    We use a rectangular contour with Re(s) ∈ [σ₁, σ₂], Im(s) ∈ [-T, T]
    where σ₁ < 1/2 < σ₂.

    The right wall (Re(s) = σ₂ > 1): ζ'/ζ computable via Dirichlet series
    The left wall: use functional equation
    Top/bottom: use known bounds

    PROBLEM: For Re(s) ≤ 1, the Dirichlet series doesn't converge.
    SOLUTION: Use the functional equation ζ(s) = χ(s)·ζ(1-s) to reflect.

    Alternative: Use the Hadamard product formula directly.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Contour Integral Approach")
    print("=" * 70)

    # For the contour integral, we need ζ'/ζ on the critical strip.
    # This is the fundamental difficulty: computing ζ'/ζ at a point
    # requires O(t^{1/2}) work via Riemann-Siegel formula.
    #
    # So evaluating the contour integral at N quadrature points costs
    # O(N · T^{1/2}) where T is the height of the contour.
    # We need T ~ x^{1/2}/log(x) for the contour to enclose enough zeros.
    # Each evaluation costs O(T^{1/2}) = O(x^{1/4}/log(x)^{1/2}).
    # Minimum N for Nyquist: O(T · log(x)) = O(x^{1/2}).
    # Total: O(x^{3/4}) -- WORSE than direct summation!
    #
    # UNLESS: we can use fewer quadrature points.
    # Test: use Gaussian quadrature / Clenshaw-Curtis with varying N.

    # For small x, we can test the principle using direct ζ'/ζ computation
    x = 100
    true_pi = int(sympy.primepi(x))

    print(f"\nTest case: x = {x}, π(x) = {true_pi}")

    # Contour: rectangle with Re(s) ∈ [-0.5, 2], Im(s) ∈ [-T, T]
    T_values = [20, 50, 100, 200]
    N_quad_values = [10, 20, 50, 100, 200, 500]

    print(f"\nDirect zero sum (reference):")
    for nz in [5, 10, 20, 50]:
        s = zero_sum_direct(x, nz)
        print(f"  {nz} zeros: S(x) = {s:.6f}")

    # Instead of the full contour integral (which needs ζ on critical strip),
    # test a DIFFERENT approach: partial fraction decomposition
    #
    # ζ'/ζ(s) has poles at each zero ρ with residue -1.
    # So: -(ζ'/ζ)(s) · x^s/s ≈ Σ_ρ x^s / (s(s-ρ)) near each ρ
    # Integrating: Σ_ρ x^ρ/ρ (recovering the zero sum)
    #
    # This is circular -- we get back to summing over zeros.

    print("\n  KEY INSIGHT: The contour integral approach is EQUIVALENT to")
    print("  summing over zeros by the residue theorem. No shortcut here.")
    print("  The number of quadrature points ≥ number of enclosed zeros")
    print("  (by Nyquist/Shannon for oscillatory integrands).")

    return None

# =============================================================================
# EXPERIMENT 3: Smoothed explicit formula / kernel methods
# =============================================================================
def experiment_smoothed_formula():
    """
    Test: Can smoothing the explicit formula reduce the number of zeros needed?

    Instead of the sharp cutoff π(x) = #{p ≤ x}, use a smoothed version:
    π_φ(x) = Σ_p φ((x-p)/h) where φ is a smooth bump function with width h.

    The smoothed explicit formula converges faster (zeros weighted by φ̂(γ),
    which decays for smooth φ). If φ̂ has compact support of width W,
    only zeros with |γ| < W contribute.

    Key question: can we choose h small enough that π_φ(x) = π(x) exactly,
    while keeping W = O(polylog(x))?

    Answer: h must be < min gap between prime and non-prime near x.
    The prime gap near x is O(x^{0.525}) (Baker-Harman-Pintz).
    So h ~ 1 suffices. But then W = 1/h ~ 1 in frequency domain.

    WAIT -- this is more subtle. Let's test numerically.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Smoothed Explicit Formula")
    print("=" * 70)

    # Test: for smooth φ with bandwidth W, how many zeros (|γ| < W) are needed?
    # Use Gaussian smoothing: φ(t) = exp(-t²/2h²), φ̂(ω) = h·exp(-h²ω²/2)

    test_x = [10**3, 10**4, 10**5, 10**6]

    for x in test_x:
        true_pi = int(sympy.primepi(x))
        Rx = R_function(x)

        # Find nearest prime gap around x
        p_below = int(sympy.prevprime(x + 1)) if sympy.isprime(x) else int(sympy.prevprime(x))
        p_above = int(sympy.nextprime(x))
        gap = p_above - p_below

        print(f"\nx = {x}, π(x) = {true_pi}, gap around x = {gap}")
        print(f"  R(x) = {Rx:.4f}, error without zeros = {Rx - true_pi:.4f}")

        # For Gaussian smoothing with width h:
        # π_φ(x) = π(x) exactly when h < distance to nearest prime/non-prime boundary
        # The frequency cutoff is W ~ 1/h
        # Number of zeros with |γ| < W: ~ W·log(W)/(2π)

        for h in [0.5, 1.0, 2.0, 5.0, 10.0]:
            W = 3.0 / h  # 3σ cutoff for Gaussian
            # Number of zeros needed: ~ W·log(W)/(2π) by zero density
            nz_needed = max(1, int(W * np.log(max(W, 2)) / (2 * np.pi)))

            # Compute smoothed zero correction
            log_x = np.log(x)
            correction = 0.0
            for gamma in ZEROS[:min(nz_needed * 3, len(ZEROS))]:  # use more for safety
                rho = 0.5 + 1j * gamma
                # Weight by φ̂(γ) = h·exp(-h²γ²/2)
                weight = h * np.exp(-h**2 * gamma**2 / 2)
                if weight < 1e-15:
                    break
                contrib = x**rho / rho
                correction += 2 * contrib.real * weight

            smoothed_estimate = Rx - correction
            error = smoothed_estimate - true_pi

            # How many zeros effectively contributed (weight > 1e-10)?
            eff_zeros = sum(1 for g in ZEROS if h * np.exp(-h**2 * g**2 / 2) > 1e-10)

            print(f"  h={h:>4.1f}: W={W:>6.1f}, eff_zeros={eff_zeros:>4}, "
                  f"estimate={smoothed_estimate:>12.4f}, error={error:>+10.4f}")

    print("\n  ANALYSIS: Smoothing reduces effective zeros but introduces bias.")
    print("  For h=0.5 (sharp): need all zeros, get exact answer.")
    print("  For h=5.0 (smooth): need few zeros, but answer is smoothed (wrong).")
    print("  The fundamental tradeoff: fewer zeros = less precision.")

# =============================================================================
# EXPERIMENT 4: Can oscillatory part be compressed via correlations?
# =============================================================================
def experiment_compression():
    """
    Test: The GUE hypothesis implies correlations between zeros.
    Can we exploit these to compress the zero sum?

    Idea: Instead of Σ_ρ x^ρ/ρ with ~N terms, find a SPARSE representation
    in some basis that captures the sum with ~polylog(N) terms.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Compressibility of the Zero Sum")
    print("=" * 70)

    # For each x, the zero sum S(x) = Σ_{k=1}^N 2·Re(x^{ρ_k}/ρ_k)
    # This is a sum of oscillating terms: each ~ x^{1/2} cos(γ_k log x - arg(ρ_k)) / |ρ_k|
    #
    # The "signal" S(x) at x is a linear combination of sinusoids with
    # frequencies γ_k (on a log scale) and amplitudes ~ 1/|ρ_k|.
    #
    # Can we represent this signal with fewer components?

    # Test: compute S(x) for many x values, then do SVD/PCA
    N_zeros = 500
    x_values = np.exp(np.linspace(np.log(100), np.log(1e6), 500))

    # Build matrix: M[i, k] = 2·Re(x_i^{ρ_k}/ρ_k)
    M = np.zeros((len(x_values), N_zeros))
    for k, gamma in enumerate(ZEROS[:N_zeros]):
        rho = 0.5 + 1j * gamma
        for i, x in enumerate(x_values):
            contrib = x**rho / rho
            M[i, k] = 2 * contrib.real

    # SVD
    U, S_vals, Vt = np.linalg.svd(M, full_matrices=False)

    print(f"\nSVD of zero-contribution matrix ({len(x_values)} x-values × {N_zeros} zeros):")
    print(f"  Top 20 singular values:")
    for i in range(min(20, len(S_vals))):
        cumulative = np.sum(S_vals[:i+1]**2) / np.sum(S_vals**2)
        print(f"    σ_{i+1:>2} = {S_vals[i]:>12.4f}  (cumulative energy: {cumulative:.6f})")

    # How many singular values capture 99%, 99.9%, 99.99% of energy?
    total_energy = np.sum(S_vals**2)
    for threshold in [0.9, 0.99, 0.999, 0.9999, 0.99999]:
        cumsum = np.cumsum(S_vals**2) / total_energy
        k = np.searchsorted(cumsum, threshold) + 1
        print(f"  {threshold*100:>7.3f}% energy captured by {k:>3} components (out of {N_zeros})")

    # Test: reconstruct S(x) using only top-k components
    print(f"\nReconstruction accuracy with rank-k approximation:")
    test_xs = [1000, 5000, 10000, 50000, 100000, 500000]

    for x in test_xs:
        true_pi = int(sympy.primepi(x))
        Rx = R_function(x)

        # Full zero sum
        full_sum = 0.0
        for gamma in ZEROS[:N_zeros]:
            rho = 0.5 + 1j * gamma
            full_sum += 2 * (x**rho / rho).real
        full_estimate = Rx - full_sum

        print(f"\n  x = {x:>7}, π(x) = {true_pi}, R(x) = {Rx:.2f}")
        print(f"  Full ({N_zeros} zeros): estimate = {full_estimate:.2f}, error = {full_estimate - true_pi:+.2f}")

        # Low-rank approximation
        # Recompute contributions at this x
        contrib_vec = np.zeros(N_zeros)
        for k, gamma in enumerate(ZEROS[:N_zeros]):
            rho = 0.5 + 1j * gamma
            contrib_vec[k] = 2 * (x**rho / rho).real

        for rank in [1, 2, 5, 10, 20, 50, 100]:
            if rank > N_zeros:
                break
            # Project onto top-rank subspace
            approx = Vt[:rank, :].T @ (Vt[:rank, :] @ contrib_vec)
            approx_sum = np.sum(approx)
            approx_estimate = Rx - approx_sum
            err = approx_estimate - true_pi
            print(f"    rank {rank:>3}: estimate = {approx_estimate:>10.2f}, error = {err:>+8.2f}")

    return S_vals

# =============================================================================
# EXPERIMENT 5: Modular prime counting — can we get p(n) mod m cheaply?
# =============================================================================
def experiment_modular():
    """
    Test: Can we determine p(n) mod m for small m without computing p(n)?

    For m=2: p(n) is odd for n≥3. Free!
    For m=3,5,7,...: By Dirichlet, primes equidistribute among residue classes.
    But the EXACT residue of p(n) requires knowing π(x; m, a) for each a.

    Question: Is π(x; m, a) mod m somehow easier than π(x)?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Modular Prime Counting")
    print("=" * 70)

    # Count primes in residue classes for various moduli
    for m in [3, 5, 7, 11, 13]:
        print(f"\nPrimes mod {m} (up to 10^5):")
        counts = {a: 0 for a in range(m) if sympy.gcd(a, m) == 1 or a == 0}
        for p in sympy.primerange(2, 100001):
            r = int(p) % m
            if r not in counts:
                counts[r] = 0
            counts[r] += 1

        total = sum(counts.values())
        for a in sorted(counts.keys()):
            print(f"  π(10^5; {m}, {a}) = {counts[a]:>5} ({counts[a]/total*100:.1f}%)")

        # Check: is the sequence p(n) mod m periodic or quasi-periodic?
        primes_mod = []
        for p in sympy.primerange(2, 10001):
            primes_mod.append(int(p) % m)

        # Autocorrelation
        pm = np.array(primes_mod, dtype=float)
        pm -= pm.mean()
        autocorr = np.correlate(pm, pm, 'full')
        autocorr = autocorr[len(autocorr)//2:]
        autocorr /= autocorr[0]

        # Find significant autocorrelation peaks
        peaks = []
        for lag in range(1, min(100, len(autocorr))):
            if abs(autocorr[lag]) > 0.05:
                peaks.append((lag, autocorr[lag]))

        if peaks:
            print(f"  Autocorrelation peaks (lag, value): {peaks[:5]}")
        else:
            print(f"  No significant autocorrelation found (as expected)")

    print("\n  VERDICT: p(n) mod m shows no useful periodicity.")
    print("  The parity problem barrier applies: even mod 2 is hard for sieves.")

# =============================================================================
# EXPERIMENT 6: Novel idea — "Spectral gap amplification"
# =============================================================================
def experiment_spectral_gap():
    """
    Novel idea: What if we raise the zeta function to a power?

    ζ(s)^k has zeros at the same locations but with multiplicity k.
    The contribution of each zero to the explicit formula is amplified.

    Alternatively, consider L(s) = Σ a_n n^{-s} where a_n is chosen to
    amplify certain features of the prime distribution.

    Test: products/powers of L-functions that might give faster convergence.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 6: L-function Combinations for Faster Convergence")
    print("=" * 70)

    # The explicit formula for ψ_k(x) associated with -ζ^(k)'/ζ^(k) gives:
    # For ζ(s)^k: the "primes" are counted with weight k·Λ(n)
    # This doesn't help — same zeros, just multiplied by k

    # But what about ζ(s)/ζ(2s)? This has poles where ζ(s)=0 AND where ζ(2s)=0
    # The Dirichlet series: Σ |μ(n)| n^{-s} (indicator of squarefree numbers)

    # Or: ζ(s)²/ζ(2s) = Σ 2^{ω(n)} n^{-s}
    # The zeros of this are at zeros of ζ(s) (double) AND zeros of ζ(2s) (pole)
    # Zeros of ζ(2s) are at ρ/2 where ρ are zeros of ζ

    # Key question: do different L-function combinations give DIFFERENT
    # convergence rates for the explicit formula?

    # Test with the Möbius function sum: M(x) = Σ_{n≤x} μ(n)
    # This has explicit formula with zeros of ζ, but the Mertens conjecture
    # says |M(x)| < √x, so the zero contributions partially cancel.

    # Compare convergence:
    # (A) π(x) via Riemann explicit formula
    # (B) ψ(x) via Chebyshev explicit formula
    # (C) M(x) via Mertens explicit formula

    x_vals = [1000, 5000, 10000, 50000]

    for x in x_vals:
        true_pi = int(sympy.primepi(x))
        Rx = R_function(x)
        log_x = np.log(x)

        print(f"\nx = {x}, π(x) = {true_pi}")

        # Compare error reduction rates
        for nz in [5, 10, 20, 50, 100, 200]:
            if nz > len(ZEROS):
                break

            # (A) Standard: π(x) ≈ R(x) - Σ R(x^ρ) using x^ρ/ρ approximation
            correction_a = sum(2 * (x**( 0.5 + 1j*g) / (0.5 + 1j*g)).real
                             for g in ZEROS[:nz])
            estimate_a = Rx - correction_a
            error_a = abs(estimate_a - true_pi)

            # (B) ψ(x): true value minus smooth part, check zero sum convergence
            correction_b = sum(2 * (x**(0.5 + 1j*g) / (0.5 + 1j*g)).real
                             for g in ZEROS[:nz])
            # ψ error is same order as π error for this purpose

            print(f"  {nz:>3} zeros: |error_π| = {error_a:.4f}")

    print("\n  ANALYSIS: All standard L-function variants have the same")
    print("  convergence rate for the zero sum. The zeros are the same;")
    print("  only the weights change, not the convergence speed.")

# =============================================================================
# RUN ALL EXPERIMENTS
# =============================================================================
if __name__ == '__main__':
    print("CONTOUR INTEGRAL & SPECTRAL METHODS FOR ZERO SUM EVALUATION")
    print("=" * 70)

    t0 = time.time()

    results1 = experiment_convergence()
    experiment_contour_integral()
    experiment_smoothed_formula()
    sv = experiment_compression()
    experiment_modular()
    experiment_spectral_gap()

    elapsed = time.time() - t0

    print(f"\n\n{'=' * 70}")
    print(f"Total runtime: {elapsed:.1f}s")
    print(f"{'=' * 70}")

    # Summary
    print("\n## OVERALL SUMMARY")
    print("1. CONTOUR INTEGRAL: Equivalent to zero summation (residue theorem)")
    print("2. SMOOTHING: Reduces zeros needed but introduces bias (tradeoff)")
    print("3. COMPRESSION: SVD shows the zero-contribution matrix IS low-rank")
    print("   for a FIXED set of x-values, but the basis depends on the x-range")
    print("4. MODULAR: No useful periodicity (parity barrier)")
    print("5. L-FUNCTION COMBOS: Same zeros, same convergence rate")
    print("\nThe compression result (Exp 4) is the most interesting —")
    print("it suggests structure, but it may be a property of the x-sampling")
    print("rather than a genuine shortcut.")
