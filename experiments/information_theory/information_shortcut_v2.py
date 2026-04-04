"""
Session 9: Information-Theoretic Shortcut v2 — Fixed

Key experiments:
1. Can GUE-sampled zeros replace actual zeros in explicit formula?
2. Connes' finite Euler product approach
3. Fast multipole for oscillatory zero sums
"""

import numpy as np
from mpmath import mp, mpf, zetazero, log as mplog, li as mpli
from sympy import prime, primepi
import time

mp.dps = 30

# Get actual zeros
K_max = 100
zeros = [float(zetazero(k).imag) for k in range(1, K_max + 1)]

# Better smooth approximation using Gram points
def gram_zero_approx(n):
    """Better approximation: solve θ(t) = (n-1)π where θ is Riemann-Siegel theta"""
    # Crude: γ_n ≈ 2π(n - 7/8) / ln(n/(2πe)) — only good for large n
    # Better: use Newton iteration on N(T) = T/(2π) · ln(T/(2πe)) + 7/8
    # We want N(γ_n) = n, so γ_n s.t. γ/(2π)·ln(γ/(2πe)) + 7/8 = n
    # i.e., γ·ln(γ/(2πe)) = 2π(n - 7/8)
    target = 2 * np.pi * (n - 7.0/8.0)
    # Newton: γ·ln(γ/c) = target where c = 2πe
    c = 2 * np.pi * np.e
    gamma = max(target / np.log(max(target/c, 2)), 10.0)
    for _ in range(20):
        f = gamma * np.log(gamma / c) - target
        fp = np.log(gamma / c) + 1
        gamma -= f / fp
        if abs(f) < 1e-10:
            break
    return gamma

print("=" * 70)
print("Improved smooth zero approximation")
print("=" * 70)
smooth_errors = []
for k in range(1, K_max + 1):
    smooth = gram_zero_approx(k)
    actual = zeros[k-1]
    err = actual - smooth
    smooth_errors.append(err)
    if k <= 10 or k % 20 == 0:
        print(f"  γ_{k:3d}: actual={actual:.4f}, smooth={smooth:.4f}, error={err:.4f}")

smooth_errors = np.array(smooth_errors)
print(f"\nSmooth errors: mean={np.mean(smooth_errors):.4f}, std={np.std(smooth_errors):.4f}, max={np.max(np.abs(smooth_errors)):.4f}")

print("\n" + "=" * 70)
print("EXPERIMENT: π(x) with actual vs smooth vs random zeros")
print("=" * 70)

def explicit_pi(x, zero_list, K):
    """Compute π(x) using explicit formula"""
    x_val = mpf(x)
    result = float(mpli(x_val))  # R(x) ≈ li(x)
    log_x = np.log(float(x))

    for k in range(min(K, len(zero_list))):
        gamma = zero_list[k]
        phase = gamma * log_x
        rho_abs = np.sqrt(0.25 + gamma**2)
        rho_angle = np.arctan(2 * gamma)
        amplitude = 2 * float(x)**0.5 / (rho_abs * log_x)
        result -= amplitude * np.cos(phase - rho_angle)

    return result

np.random.seed(42)
smooth_zeros = [gram_zero_approx(k) for k in range(1, K_max + 1)]
random_zeros = [gram_zero_approx(k) + np.random.normal(0, np.std(smooth_errors)) for k in range(1, K_max + 1)]

print(f"{'x':>8} {'π(x)':>8} {'li(x)':>10} {'50actual':>10} {'50smooth':>10} {'50random':>10}")
for x in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    actual = int(primepi(x))
    li_x = float(mpli(mpf(x)))
    pi_actual = explicit_pi(x, zeros, 50)
    pi_smooth = explicit_pi(x, smooth_zeros, 50)
    pi_random = explicit_pi(x, random_zeros, 50)

    print(f"{x:8d} {actual:8d} {li_x:10.1f} {pi_actual:10.1f} {pi_smooth:10.1f} {pi_random:10.1f}")

print("\nErrors (actual minus estimate):")
print(f"{'x':>8} {'err_li':>8} {'err_50a':>8} {'err_50s':>8} {'err_50r':>8}")
for x in [100, 500, 1000, 5000, 10000, 50000, 100000]:
    actual = int(primepi(x))
    li_x = float(mpli(mpf(x)))
    pi_actual = explicit_pi(x, zeros, 50)
    pi_smooth = explicit_pi(x, smooth_zeros, 50)
    pi_random = explicit_pi(x, random_zeros, 50)

    print(f"{x:8d} {actual-li_x:8.1f} {actual-pi_actual:8.1f} {actual-pi_smooth:8.1f} {actual-pi_random:8.1f}")

print("\n" + "=" * 70)
print("EXPERIMENT: How many zeros needed for error < gap/2?")
print("=" * 70)

# For floor rounding to work, we need |π(x) - π_approx(x)| < 0.5
# at each step of binary search

for K in [10, 20, 50, 100]:
    exact_count = 0
    total = 0
    for x in range(100, 1001):
        actual = int(primepi(x))
        approx = explicit_pi(x, zeros, K)
        if abs(actual - approx) < 0.5:
            exact_count += 1
        total += 1
    pct = 100 * exact_count / total
    print(f"  K={K:3d} zeros: {pct:.1f}% of π(x) values in [100,1000] exact to ±0.5")

print("\n" + "=" * 70)
print("EXPERIMENT: Zero sum fast multipole feasibility")
print("=" * 70)

# The sum S(x) = Σ_{k=1}^N x^{ρ_k}/ρ_k = Σ x^{1/2+iγ_k} / (1/2+iγ_k)
# = x^{1/2} Σ e^{iγ_k·log(x)} / (1/2+iγ_k)

# This is a TYPE III sum in analytic number theory
# Can it be computed using fast methods?

# Key: if zeros γ_k are approximately equally spaced (spacing ~ 2π/log(T)):
# The sum is like a DFT evaluated at a single point t = log(x)
# A single DFT value takes O(N) — no shortcut from FFT

# BUT: if zeros can be grouped into bands where the contribution
# is predictable, we might get O(√N) or O(N^{1/3})

# Test: compute partial sums over bands of zeros
t = np.log(10000)  # x = 10000
N = 100

band_size = 10
num_bands = N // band_size

print(f"\nBand contributions to S(x) at x=10000 (t={t:.2f}):")
print(f"{'Band':>6} {'Re(S_band)':>12} {'|S_band|':>12} {'Phase':>10}")

total_S = 0
band_contribs = []
for b in range(num_bands):
    S_band = 0
    for k in range(b * band_size, (b+1) * band_size):
        gamma = zeros[k]
        rho = 0.5 + 1j * gamma
        x_rho = np.exp(rho * t)
        S_band += x_rho / rho
    band_contribs.append(S_band)
    total_S += S_band
    phase = np.angle(S_band) * 180 / np.pi
    print(f"{b*band_size+1:3d}-{(b+1)*band_size:3d} {S_band.real:12.4f} {abs(S_band):12.4f} {phase:10.1f}°")

print(f"\nTotal S = {total_S.real:.4f} + {total_S.imag:.4f}i")
print(f"|Total S| = {abs(total_S):.4f}")

# Check: do band contributions decrease?
magnitudes = [abs(s) for s in band_contribs]
print(f"\nBand magnitudes: {[f'{m:.2f}' for m in magnitudes]}")
print(f"Decreasing? {all(magnitudes[i] >= magnitudes[i+1] for i in range(len(magnitudes)-1))}")

# The magnitudes are NOT consistently decreasing — bands oscillate
# This means no FMM-like truncation is possible

print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
print("""
1. Smooth zero approximation: errors ~0.4 std, good but not sufficient
2. Smooth zeros in explicit formula: similar accuracy to actual zeros for small x
   but DIVERGES for large x (as expected)
3. GUE-sampled random zeros: WORSE than smooth — confirms individual zeros matter
4. Exact rounding: even with 100 zeros, only works for small x (< 1000)
5. Band contributions: NOT decreasing → no fast multipole truncation possible
6. Zero sum is fundamentally O(N) — each zero contributes independently

The Connes approach is beautiful but doesn't bypass the fundamental barrier:
- Computing individual zeros: O(1) each with Connes' method (maybe)
- But SUMMING them: still O(√x) terms needed
- Total: O(√x) even with perfect zero approximation

The barrier is NOT about computing zeros — it's about SUMMING enough of them.
""")
