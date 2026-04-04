"""
Session 9: Information-Theoretic Shortcut

RADICAL NEW IDEA: "Compressed Explicit Formula"

The explicit formula π(x) = R(x) - Σ_ρ R(x^ρ) needs O(√x) zeros.
But what if we don't need INDIVIDUAL zeros — just AGGREGATE statistics?

Key insight from random matrix theory:
- Zeta zeros are GUE-distributed locally
- The n-point correlations are KNOWN (Montgomery, Hejhal, Rudnick-Sarnak)
- The pair correlation is: 1 - (sin(πu)/(πu))² (Montgomery)

For the explicit formula sum S(x) = Σ_ρ R(x^ρ):
- Each term ≈ x^{1/2} · e^{iγ·log(x)} / (ρ · log(x))
- This is like Σ e^{iγ_k · t} where t = log(x)

By Weyl's equidistribution + GUE statistics:
- The sum S(x) = Σ e^{iγ_k·t} for t = log(x)
- This is the SPECTRAL FORM FACTOR of the GUE ensemble
- At time t: |S(t)|² ~ t for t < t_H = 2π·N(T)/T (Heisenberg time)
- After t_H: |S(t)|² saturates to N(T)

The Heisenberg time for x zeros up to height T is t_H ~ T.
For x ~ 10^102: we need T ~ √x ~ 10^51, so t_H ~ 10^51.
log(x) ~ 235 << 10^51, so we're in the early-time regime.

In the early-time regime: S(t) ~ √(t/2π) · random_complex
This gives |π(x) - R(x)| ~ √(log(x)/2π) · √x ≈ √(235/6.28) · 10^51
This is ~6.1 × 10^51 — consistent with known error bounds.

BUT: the PHASE of S(t) is what matters for the exact value of π(x).
The phase is uniformly distributed → no shortcut.

Let me verify this numerically and explore if there's a way around it.
"""

import numpy as np
from mpmath import mp, mpf, zetazero, log, exp, pi, li, sqrt
from sympy import prime, primepi
import time

mp.dps = 30

print("=" * 70)
print("EXPERIMENT 1: Spectral Form Factor of Zeta Zeros")
print("=" * 70)

# Compute the sum S(t) = Σ_{k=1}^K e^{iγ_k·t} for various t
K_max = 100  # Number of zeros to use
print(f"Computing first {K_max} zeta zeros...")
zeros = []
for k in range(1, K_max + 1):
    gamma = float(zetazero(k).imag) if hasattr(zetazero(k), 'imag') else float(zetazero(k))
    zeros.append(gamma)
print(f"Done. γ_1 = {zeros[0]:.4f}, γ_{K_max} = {zeros[-1]:.4f}")

# Form factor at various "times" t = log(x) for different x
print("\nSpectral form factor |S(t)|² vs t:")
for x_exp in [2, 3, 4, 5, 6, 8, 10, 15, 20]:
    x = 10.0 ** x_exp
    t = np.log(x)
    S = sum(np.exp(1j * g * t) for g in zeros)
    form_factor = abs(S)**2 / K_max  # Normalized
    actual_pi = primepi(int(x)) if x_exp <= 8 else None
    R_x = float(li(mpf(x)))  # Rough R(x) ≈ li(x)
    if actual_pi is not None:
        error = actual_pi - R_x
        print(f"  x=10^{x_exp:2d}: t={t:.1f}, |S|²/K={form_factor:.2f}, π(x)={actual_pi}, R(x)≈{R_x:.0f}, error={error:.0f}")
    else:
        print(f"  x=10^{x_exp:2d}: t={t:.1f}, |S|²/K={form_factor:.2f}")

print("\n" + "=" * 70)
print("EXPERIMENT 2: Can GUE Statistics Replace Individual Zeros?")
print("=" * 70)

# The GUE prediction for the pair correlation function is:
# R_2(u) = 1 - (sin(πu)/(πu))²

# The explicit formula sum can be rewritten using pair correlations:
# Var[S(t)] = K + 2·Σ_{j<k} cos((γ_j - γ_k)·t)
# ≈ K + 2·K²·∫ R_2(u)·cos(2π·u·t/Δ)·du where Δ is mean spacing

# Mean spacing of zeros near height T: Δ ≈ 2π/log(T/(2π))
mean_spacing = 2 * np.pi / np.log(zeros[-1] / (2 * np.pi))
print(f"Mean spacing of zeros near γ={zeros[-1]:.1f}: Δ ≈ {mean_spacing:.4f}")

# Actual pair differences
diffs = []
for i in range(len(zeros)):
    for j in range(i+1, min(i+20, len(zeros))):
        diffs.append((zeros[j] - zeros[i]) / mean_spacing)

diffs = np.array(diffs)
print(f"Number of pair differences computed: {len(diffs)}")

# Compare with GUE prediction
bins = np.linspace(0, 5, 50)
hist, _ = np.histogram(diffs, bins=bins, density=True)
bin_centers = (bins[:-1] + bins[1:]) / 2
gue_pred = 1 - (np.sin(np.pi * bin_centers) / (np.pi * bin_centers + 1e-10))**2

# Compute L1 error
l1_error = np.sum(np.abs(hist - gue_pred) * (bins[1] - bins[0]))
print(f"L1 distance from GUE pair correlation: {l1_error:.4f}")

print("\n" + "=" * 70)
print("EXPERIMENT 3: Approximating π(x) Using GUE-Sampled Zeros")
print("=" * 70)

# IDEA: Instead of actual zeros, use RANDOM zeros with GUE statistics.
# If this gives the same accuracy as actual zeros, then the individual
# zero values don't matter — only their statistics do.

# Generate GUE-like zeros using the smooth approximation + random perturbation
# γ_n ≈ 2πn/log(n/(2πe)) (Backlund smooth approximation)
# Add GUE-distributed perturbation

def smooth_zero(n):
    """Backlund smooth approximation for nth zeta zero"""
    if n <= 0:
        return 0
    # Improved: γ_n ≈ 2πn / W(n/e) where W is Lambert W
    # Simpler: γ_n ≈ 2π(n - 7/8) / log((n - 7/8)/(2πe))
    n_adj = n - 7.0/8.0
    if n_adj <= 0:
        n_adj = 0.1
    return 2 * np.pi * n_adj / np.log(n_adj / (2 * np.pi * np.e))

# Compare smooth approximation with actual zeros
print("Smooth approximation vs actual zeros:")
smooth_errors = []
for k in range(1, K_max + 1):
    smooth = smooth_zero(k)
    actual = zeros[k-1]
    err = actual - smooth
    smooth_errors.append(err)
    if k <= 10 or k % 20 == 0:
        print(f"  γ_{k:3d}: actual={actual:.4f}, smooth={smooth:.4f}, error={err:.4f}")

smooth_errors = np.array(smooth_errors)
print(f"\nSmooth approximation errors:")
print(f"  Mean: {np.mean(smooth_errors):.4f}")
print(f"  Std: {np.std(smooth_errors):.4f}")
print(f"  Max |error|: {np.max(np.abs(smooth_errors)):.4f}")

# Now: compute π(x) using actual zeros vs smooth zeros vs GUE-sampled zeros
def explicit_pi_from_zeros(x, zero_list, K):
    """Compute π(x) approximation using explicit formula with given zeros"""
    x = float(x)
    log_x = np.log(x)
    result = float(li(mpf(x)))  # R(x) ≈ li(x)

    for k in range(min(K, len(zero_list))):
        gamma = zero_list[k]
        # Contribution of pair ρ, ρ̄:
        # -2·Re(li(x^{1/2+iγ})) ≈ -2·x^{1/2}·cos(γ·log(x)) / ((1/4+γ²)^{1/2}·log(x))
        phase = gamma * log_x
        # More accurate: -2·Re(Ei((1/2+iγ)·log(x))) / log(x)
        # Simplified: -2·x^{0.5}·cos(phase) / (gamma * log_x) approximately
        amplitude = 2 * x**0.5 / (np.sqrt(0.25 + gamma**2) * log_x)
        result -= amplitude * np.cos(phase - np.arctan(2*gamma))

    return result

print("\n\nπ(x) using actual vs smooth vs random zeros:")
print(f"{'x':>10} {'actual π(x)':>12} {'K actual':>12} {'K smooth':>12} {'K random':>12}")

np.random.seed(42)
random_zeros = [smooth_zero(k) + np.random.normal(0, 0.5) for k in range(1, K_max + 1)]

for x in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    actual = primepi(x)
    K = 50
    pi_actual_zeros = explicit_pi_from_zeros(x, zeros, K)
    pi_smooth_zeros = explicit_pi_from_zeros(x, [smooth_zero(k) for k in range(1, K_max+1)], K)
    pi_random_zeros = explicit_pi_from_zeros(x, random_zeros, K)

    print(f"{x:10d} {actual:12d} {pi_actual_zeros:12.1f} {pi_smooth_zeros:12.1f} {pi_random_zeros:12.1f}")

print("\n" + "=" * 70)
print("EXPERIMENT 4: Connes-Inspired Finite Euler Product")
print("=" * 70)

# Connes' key idea: use the Weil explicit formula with FINITE Euler product
# over small primes to constrain the zero locations.
# The Weil explicit formula for a test function f:
# Σ_ρ f̂(γ_ρ) = f(0)·(log 4π + γ) - ∫... + Σ_{p≤P} Σ_k f̂(k·log p)/p^{k/2}

# With P=13 (primes 2,3,5,7,11,13), the RHS is computable.
# The LHS constrains the zero locations.

# Let's try: for a test function f_t(x) = e^{-x²/(2σ²)} (Gaussian),
# f̂_t(ω) = σ√(2π) · e^{-σ²ω²/2}
# The sum Σ_ρ f̂(γ_ρ) = σ√(2π) · Σ_ρ e^{-σ²γ_ρ²/2}

# For large σ: this samples many zeros equally
# For small σ: this focuses on low-lying zeros

# The prime sum: Σ_{p≤13} Σ_k σ√(2π)·e^{-σ²(k·log p)²/2} / p^{k/2}

small_primes = [2, 3, 5, 7, 11, 13]

def weil_prime_sum(sigma, P_list=small_primes, k_max=10):
    """Compute the prime-side of Weil explicit formula with Gaussian test function"""
    total = 0
    for p in P_list:
        for k in range(1, k_max + 1):
            log_p = np.log(p)
            total += sigma * np.sqrt(2*np.pi) * np.exp(-sigma**2 * (k*log_p)**2 / 2) / p**(k/2)
    return total

def weil_zero_sum(sigma, zero_list):
    """Compute the zero-side of Weil explicit formula with Gaussian test function"""
    total = 0
    for gamma in zero_list:
        total += sigma * np.sqrt(2*np.pi) * np.exp(-sigma**2 * gamma**2 / 2)
    total *= 2  # Both ρ and ρ̄
    return total

def weil_constant_term(sigma):
    """The constant/continuous part of Weil formula"""
    gamma_euler = 0.5772156649
    return sigma * np.sqrt(2*np.pi) * (np.log(4*np.pi) + gamma_euler)

print("Weil formula verification (Gaussian test function):")
for sigma in [0.05, 0.1, 0.2, 0.5, 1.0]:
    prime_side = weil_prime_sum(sigma)
    zero_side = weil_zero_sum(sigma, zeros[:50])
    constant = weil_constant_term(sigma)
    # The integral term is more complex, skip for now
    rhs = prime_side + constant
    print(f"  σ={sigma:.2f}: prime_sum={prime_side:.4f}, zero_sum={zero_side:.4f}, const={constant:.4f}, RHS≈{rhs:.4f}")

# The key question: can we INVERT this to get zeros from primes?
# For each σ, the equation Σ_ρ f̂(γ_ρ) = known_RHS gives one constraint.
# With many σ values, we get many constraints → can recover zeros.

# This is essentially a MOMENT problem: from moments, recover the measure.
# The measure is the spectral measure Σ_ρ δ(γ - γ_ρ).
# The moments are: M_k = Σ_ρ γ_ρ^k (power sums of zeros).

# Can we get moments from the Weil formula?
# Using f(x) = x^{2k} · e^{-x²}: f̂(ω) involves Hermite polynomials
# Each such f gives Σ_ρ H_{2k}(γ_ρ) · e^{-γ_ρ²} = known(primes)

# From these moments/Hermite sums, we can reconstruct zeros.
# This is the Connes approach!

# Let's try it numerically with a simplified version:
print("\n\nMoment recovery of zeros from small primes:")
print("Computing moments of γ_ρ from Weil formula...")

# Use Li's criterion: λ_n = Σ_ρ [1 - (1-1/ρ)^n]
# RH iff λ_n > 0 for all n ≥ 1
# These are related to power sums of zeros via Stirling numbers

# Simpler: direct power sums S_k = Σ_ρ 1/ρ^k
# Known from the Laurent expansion of ξ'/ξ(s)

# The key identity: Σ_ρ 1/ρ^k = (coefficient in expansion of -ζ'/ζ)
# For k=1: Σ 1/ρ = (γ/2 + 1 - log(4π)/2 + ...)

# Actually, let's compute Σ γ^{2k} for our 100 zeros and see if it matches
# any predictable formula
for k in [1, 2, 3, 4]:
    moment = sum(g**(2*k) for g in zeros)
    print(f"  Σ γ^{2*k} (first 100 zeros) = {moment:.2f}")
    # These grow very fast — no simple pattern

print("\n" + "=" * 70)
print("EXPERIMENT 5: Fast π(x) via Connes' Finite Approximation")
print("=" * 70)

# The practical question: if we approximate zeros using only small primes
# (Connes' method), how accurate is π(x)?

# Step 1: Use the smooth approximation for zeros (which uses NO primes)
# Step 2: Correct using the Weil formula constraints from small primes
# Step 3: Use corrected zeros in the explicit formula

# Step 1 is done above. For step 2:
# The Weil formula gives us constraints like:
# Σ_ρ e^{-σ²γ²/2} = F(σ, primes≤13)
# We can use these to refine the smooth approximation

# Let's measure: does adding Weil constraints improve π(x)?

# Baseline: π(x) using only R(x) (no zeros)
# Level 1: π(x) using smooth zeros
# Level 2: π(x) using actual first K zeros

print(f"{'x':>8} {'π(x)':>8} {'R(x)':>10} {'err_R':>8} {'smooth_50':>10} {'err_s':>8} {'actual_50':>10} {'err_a':>8}")
for x in [100, 500, 1000, 5000, 10000, 50000]:
    actual = primepi(x)
    r_x = float(li(mpf(x)))
    smooth_list = [smooth_zero(k) for k in range(1, 51)]
    pi_smooth = explicit_pi_from_zeros(x, smooth_list, 50)
    pi_actual = explicit_pi_from_zeros(x, zeros, 50)

    print(f"{x:8d} {actual:8d} {r_x:10.1f} {r_x-actual:8.1f} {pi_smooth:10.1f} {pi_smooth-actual:8.1f} {pi_actual:10.1f} {pi_actual-actual:8.1f}")

print("\n" + "=" * 70)
print("CRITICAL ANALYSIS")
print("=" * 70)
print("""
The Connes approach is mathematically beautiful but faces the same barrier:

1. The Weil formula with finite primes gives CONSTRAINTS on zeros,
   but not enough to DETERMINE them for use in π(x) for large x.

2. For π(x) with x ~ 10^102, we need the sum Σ R(x^ρ) over ~10^51 zeros.
   Even if each zero is approximated to 10^{-55} accuracy (Connes' claim),
   the cumulative error from 10^51 terms is 10^51 × 10^{-55} = 10^{-4}.
   This might actually be small enough! But we need 10^51 zeros computed.

3. The bottleneck SHIFTS: from "compute zeros" to "sum 10^51 terms".
   Even if each term is O(1) to compute, summing 10^51 terms takes 10^51 ops.

4. Unless: the sum can be computed in batch using the STATISTICAL PROPERTIES
   of the approximated zeros. This would require:
   - Fast multipole method (FMM) for the oscillatory sum
   - Or: the sum telescopes due to GUE statistics
   - Neither is known to work.

CONCLUSION: Connes' method potentially helps with INDIVIDUAL zero accuracy
but does NOT reduce the NUMBER of zeros needed. The barrier remains:
O(√x) zeros × O(1) per zero = O(√x) total.
""")
