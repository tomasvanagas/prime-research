#!/usr/bin/env python3
"""
Session 10: Random Matrix Theory + Zeta Zero Proxy

IDEA: Use GUE random matrix eigenvalues as proxies for zeta zeros
to compute the correction term Σ_ρ R(x^ρ) statistically.

Also: test whether the correction term can be bounded tightly enough
that we can identify p(n) EXACTLY from the bounded window.
"""

import numpy as np
import math
from sympy import prime, primepi
from mpmath import mp, mpf, li, log as mplog, pi as mpi, zeta, zetazero

mp.dps = 30

# ============================================================
# PART 1: Compute actual zeta zeros and their contribution
# ============================================================
print("=" * 60)
print("PART 1: Actual Zeta Zero Contributions")
print("=" * 60)

# Compute first K zeta zeros
K = 50
print(f"Computing first {K} zeta zeros...")
zeros = []
for k in range(1, K + 1):
    z = zetazero(k)
    zeros.append(float(z.imag))

print(f"First 10 zeros (imaginary parts): {[f'{g:.4f}' for g in zeros[:10]]}")

# For each test n, compute:
# 1. The exact p(n)
# 2. R^{-1}(n) approximation
# 3. Contribution of first K zeros

def R_function(x):
    """Riemann R function"""
    result = mpf(0)
    for k in range(1, 100):
        mu_k = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1][k] if k <= 10 else 0
        if mu_k == 0 and k > 10:
            continue
        if k <= 10 and mu_k != 0:
            result += mpf(mu_k) / k * li(x ** (mpf(1)/k))
    return result

def R_inv(n):
    """Inverse R function"""
    x = mpf(n) * mplog(mpf(n))
    for _ in range(50):
        rx = li(x)  # Approximate R(x) ≈ li(x) for Newton step
        if abs(rx - n) < mpf('1e-20'):
            break
        x = x - (rx - n) * mplog(x)
    return x

# Compute zero contributions for test values
print("\nZero contributions to π(x):")
for n_test in [1000, 5000, 10000]:
    p_n = prime(n_test)
    x = float(R_inv(n_test))

    # Σ_ρ R(x^ρ) ≈ Σ_k 2*Re(li(x^{1/2+iγ_k}))/(number_of_Mobius_terms)
    # Simplified: contribution ≈ -Σ_k 2*cos(γ_k * log(x)) / (γ_k * sqrt(x) * log(x))
    # Actually: R(x^ρ) ≈ li(x^ρ) for the dominant term
    # li(x^{1/2+iγ}) ≈ x^{1/2} * e^{iγ*log(x)} / ((1/2+iγ)*log(x))

    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    total_contribution = 0
    for gamma in zeros:
        # Real part of li(x^{1/2+iγ})
        phase = gamma * log_x
        denom_real = 0.5
        denom_imag = gamma
        denom_sq = denom_real**2 + denom_imag**2

        # x^{1/2+iγ} = sqrt(x) * (cos(γ*log(x)) + i*sin(γ*log(x)))
        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)

        # li(x^s) ≈ x^s / (s * log(x)) for large x
        # Real part: sqrt(x)/(log(x)*|s|^2) * (s_real*cos - s_imag*sin)
        contrib = sqrt_x / (log_x * denom_sq) * (denom_real * cos_phase + denom_imag * sin_phase)
        total_contribution += 2 * contrib  # Factor of 2 for conjugate pair

    actual_error = p_n - round(x)
    print(f"  n={n_test}: p(n)={p_n}, R^-1≈{x:.0f}, "
          f"Σ_{K}_zeros≈{total_contribution:.1f}, actual_error={actual_error}")

# ============================================================
# PART 2: GUE Random Matrix Proxy
# ============================================================
print("\n" + "=" * 60)
print("PART 2: GUE Random Matrix Proxy for Zeta Zeros")
print("=" * 60)

# The N×N GUE has eigenvalues with spacing ~π/N at the origin
# Zeta zeros have spacing ~2π/log(T/(2π)) near height T
# To match T zeros up to height T, use N ≈ T*log(T)/(2π)

# Generate GUE eigenvalues
def gue_eigenvalues(N):
    """Generate eigenvalues of N×N GUE matrix"""
    # GUE: H = (A + A†)/2 where A has complex Gaussian entries
    A = (np.random.randn(N, N) + 1j * np.random.randn(N, N)) / np.sqrt(2*N)
    H = (A + A.conj().T) / 2
    eigenvalues = np.linalg.eigvalsh(H)
    return eigenvalues

# Compare GUE spacing statistics with zeta zero spacing
zero_spacings = np.diff(zeros)
mean_spacing = np.mean(zero_spacings)
normalized_spacings = zero_spacings / mean_spacing

print(f"Zeta zero spacings (normalized):")
print(f"  Mean: {np.mean(normalized_spacings):.4f}")
print(f"  Std:  {np.std(normalized_spacings):.4f}")
print(f"  Min:  {np.min(normalized_spacings):.4f}")

# Generate GUE and compare
N_gue = 100
gue_eigs = gue_eigenvalues(N_gue)
# Take central eigenvalues (bulk of spectrum)
central = gue_eigs[N_gue//4:3*N_gue//4]
gue_spacings = np.diff(central)
gue_mean_spacing = np.mean(gue_spacings)
gue_normalized = gue_spacings / gue_mean_spacing

print(f"\nGUE ({N_gue}×{N_gue}) spacing distribution:")
print(f"  Mean: {np.mean(gue_normalized):.4f}")
print(f"  Std:  {np.std(gue_normalized):.4f}")

# The Wigner surmise: P(s) = (π/2)s * exp(-πs²/4)
# Both should follow this approximately

# KEY TEST: If we use GUE eigenvalues AS IF they were zeta zeros,
# how well does the resulting π(x) approximation work?

def pi_from_zeros(x, gamma_list, K_use):
    """Approximate π(x) using explicit formula with given 'zeros'"""
    # π(x) ≈ R(x) - Σ_k R(x^{ρ_k})
    # ≈ li(x) - Σ_k 2*Re(li(x^{1/2+iγ_k}))

    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # Main term: li(x)
    result = float(li(mpf(x)))

    # Zero contributions
    for gamma in gamma_list[:K_use]:
        phase = gamma * log_x
        s_real, s_imag = 0.5, gamma
        s_sq = s_real**2 + s_imag**2

        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)

        contrib = sqrt_x / (log_x * s_sq) * (s_real * cos_phase + s_imag * sin_phase)
        result -= 2 * contrib

    return result

# Test: actual zeros vs GUE proxy
print("\nπ(x) estimates using actual zeros vs GUE proxy:")
# Scale GUE eigenvalues to match zeta zero density
# Zeta zeros near height T have density log(T/(2π))/(2π)
# First K=50 zeros go up to γ_50 ≈ 143
T_max = zeros[-1]
gue_scaled = np.sort(np.abs(gue_eigs[:K]))  # Take K eigenvalues
# Rescale to [0, T_max]
if len(gue_scaled) > 0:
    gue_scaled = list(gue_scaled * T_max / gue_scaled.max())

for x_test in [1000, 10000, 100000]:
    pi_actual = int(primepi(x_test))
    pi_real_zeros = pi_from_zeros(x_test, zeros, K)
    pi_gue_zeros = pi_from_zeros(x_test, list(gue_scaled), min(K, len(gue_scaled)))
    pi_no_zeros = float(li(mpf(x_test)))

    print(f"  x={x_test}: actual π={pi_actual}, "
          f"li(x)={pi_no_zeros:.1f}, "
          f"real zeros={pi_real_zeros:.1f}, "
          f"GUE proxy={pi_gue_zeros:.1f}")

# ============================================================
# PART 3: Statistical bounds on the correction
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Statistical Bounds on Correction")
print("=" * 60)

# Under RH: |π(x) - li(x)| ≤ (1/(8π)) √x log(x) for x ≥ 2657
# This gives an explicit bound on the error

# For x = p(10^100) ≈ 2.3 × 10^102:
# Bound = (1/(8π)) * 10^51 * 236 ≈ 10^52
# The actual error is WITHIN this bound

# But can we NARROW the bound using partial zero information?
# If we know the first K zeros exactly, the remaining error is:
# |Σ_{k>K} R(x^{ρ_k})| ≤ ???

# Estimate: each term contributes O(x^{1/2}/γ_k)
# So remaining error ≈ x^{1/2} * Σ_{k>K} 1/γ_k
# ≈ x^{1/2} * ∫_K^∞ 1/t · (log(t/(2π))/(2π)) dt  (using zero density)
# This integral diverges logarithmically!

# But with CANCELLATION: the terms oscillate, so by Dirichlet test:
# |Σ_{k>K} cos(γ_k * log(x)) / γ_k| ≈ 1/γ_K * (log(x) factor)

# So the error after K zeros is approximately:
# E(K) ≈ x^{1/2} / (γ_K * log(x))

for K_test in [10, 50, 100, 1000, 10000]:
    gamma_K = K_test * 2 * math.pi / math.log(K_test + 10)  # Approximate γ_K

    # For x = 10^102 (p(10^100)):
    x_exp = 102  # log10(x)
    sqrt_x = 10**(x_exp/2)  # = 10^51
    log_x = x_exp * math.log(10)  # ≈ 235

    error_bound = sqrt_x / (gamma_K * log_x)
    bits_remaining = math.log2(error_bound) if error_bound > 1 else 0

    print(f"  K={K_test:>5d}: γ_K≈{gamma_K:.0f}, "
          f"error bound ≈ 10^{math.log10(max(error_bound,1)):.1f}, "
          f"bits remaining ≈ {bits_remaining:.0f}")

print("\nTo reduce error below 1 (exact answer), need K such that:")
print("  x^{1/2} / (γ_K * log(x)) < 1")
print("  γ_K > x^{1/2} / log(x) ≈ 10^51 / 235 ≈ 4×10^48")
print("  K ≈ γ_K * log(γ_K) / (2π) ≈ 10^49")
print("  This is the SUMMATION BARRIER: need ~10^49 zeros")

# ============================================================
# PART 4: Can Monte Carlo sampling of zeros work?
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Monte Carlo Sampling of Zero Sum")
print("=" * 60)

# Instead of summing ALL ~10^49 zeros, can we SAMPLE?
# If we randomly sample M zeros from the pool of ~10^49:
# E[sum] = (10^49/M) * Σ_{sampled} R(x^{ρ_k})
# Var[sum] = ???

# The terms R(x^{ρ_k}) have magnitude ~x^{1/2}/γ_k and oscillate
# For a sum of N random terms of magnitude σ, the std is σ√N
#
# Individual term magnitude: ~x^{1/2}/γ_k ≈ 10^51/γ_k
# For γ_k uniform in [0, 10^49]: average magnitude ≈ 10^51/10^49 = 100
#
# Sum of 10^49 such terms: mean ≈ actual value, std ≈ 100 * √(10^49) = 10^{26.5}
# Actual sum ≈ 10^{51} (one-sigma from zero)
#
# If we sample M terms and scale:
# Estimated sum = (10^49/M) * Σ_{sampled}
# Std of estimate = 100 * √(10^49/M) * (10^49/M)
#                 = 100 * 10^{49/2} * 10^{49} / M^{3/2}
#
# For this to be < 1: M^{3/2} > 10^{75.5}, so M > 10^{50.3}
# Still astronomical!

print("Monte Carlo sampling analysis:")
print("  To estimate Σ_{k=1}^{10^49} R(x^{ρ_k}) with error < 1:")
print("  Need ~10^50 samples — NO improvement over direct summation")
print("  Reason: high variance from oscillatory terms prevents sampling")

# But what about importance sampling? Weight by |R(x^{ρ_k})|?
# The largest contributions come from small γ_k.
# First ~√x zeros contribute O(1) each.
# Last zeros contribute O(x^{1/2}/γ_K) < O(1) each.
#
# If we sample more heavily from small γ, we improve variance...
# But the number of "important" zeros is still ~x^{1/2}/log(x) ≈ 10^49

print("\n  Importance sampling:")
print("  'Important' zeros (contributing O(1) each): first ~10^49")
print("  ALL of them are important — no subset suffices")

# ============================================================
# PART 5: Novel idea - can we compute the SUM without individual terms?
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Can We Compute the Sum Without Individual Terms?")
print("=" * 60)

# The explicit formula says:
# π(x) = li(x) - Σ_ρ li(x^ρ) + ∫_x^∞ dt/(t(t²-1)log(t)) - log(2)
#
# The sum Σ_ρ li(x^ρ) can be written as a contour integral:
# Σ_ρ li(x^ρ) = -(1/2πi) ∫_{c-i∞}^{c+i∞} (ζ'/ζ)(s) * li(x^s) ds
#
# This integral passes through the pole at s=1 (residue li(x))
# and captures all the zeros.
#
# Can we evaluate this integral by SADDLE POINT / STEEPEST DESCENT?

# The integrand is: F(s) = -(ζ'/ζ)(s) * x^s / (s * log(x))
# For x large, x^s is HUGE for Re(s) > 0 and oscillates wildly
# The "saddle point" would be where d/ds[s*log(x) + log((ζ'/ζ)(s))] = 0
# This is: log(x) + (ζ''/ζ' - (ζ'/ζ)')(s) = 0
# For s on the critical line: the second term oscillates rapidly
# No useful saddle point exists for x > 10^10 (too oscillatory)

print("Contour integral for Σ_ρ li(x^ρ):")
print("  F(s) = -(ζ'/ζ)(s) * x^s / (s·log(x))")
print("  Saddle point: log(x) = -(ζ''/ζ')(s) + ((ζ'/ζ)')²(s)")
print("  For x = 10^100: log(x) = 230, but ζ'/ζ oscillates with amplitude ~log(T)")
print("  NO useful saddle point — integration requires resolving ALL oscillations")
print("  This is equivalent to computing ALL zeros → back to summation barrier")

# ============================================================
print("\n" + "=" * 60)
print("OVERALL CONCLUSION")
print("=" * 60)
print("""
Random Matrix Theory / Statistical approaches:
1. GUE proxy gives WRONG specific values (right statistics, wrong phases)
2. Monte Carlo on zeros needs ~10^50 samples (no improvement)
3. Importance sampling doesn't help (all 10^49 zeros are "important")
4. Saddle-point integration of the explicit formula fails (too oscillatory)
5. The contour integral IS the explicit formula — no shortcut

The correction term Σ_ρ R(x^ρ) is:
- Deterministic (not random)
- Non-compressible (spectral flatness ~0.9)
- Requires ~10^49 independent pieces of information
- Cannot be sampled, averaged, or statistically estimated

This is the INFORMATION-THEORETIC CORE of the barrier.
""")
