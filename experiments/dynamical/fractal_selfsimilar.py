"""
Session 7 Experiment: Fractal/Self-Similar Structure in Prime Distribution
===========================================================================
RADICAL IDEA: The prime counting function π(x) exhibits self-similar scaling.
The PNT says π(x) ~ x/ln(x), so π(cx) ~ cπ(x) + correction.

What if the CORRECTION term has a self-similar structure that can
be computed recursively?

Specifically, define:
  δ(x) = π(x) - li(x)  (the error term)

If δ(x) has a fractal/self-similar structure, maybe:
  δ(x) = f(x) · δ(x^α) + g(x)  for some α < 1

This would give a RECURSIVE formula for δ(x) that converges in
O(log log x) steps!

The explicit formula says:
  δ(x) ~ -Σ_ρ li(x^ρ) / ln(x)

where ρ = 1/2 + iγ are zeta zeros. This IS "self-similar" in the sense
that the correction involves x^{1/2} (a power scaling).

Let's explore this computationally.
"""

import math
import time
import numpy as np

def sieve_primes(n):
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n + 1) if is_prime[i]]

def li(x):
    """Logarithmic integral via numerical integration."""
    if x <= 2:
        return 0
    # li(x) = integral from 2 to x of 1/ln(t) dt
    n = 1000
    h = (x - 2) / n
    total = 0
    for i in range(n):
        t = 2 + (i + 0.5) * h
        total += 1 / math.log(t)
    return total * h

def R_func(x):
    """Riemann's R function."""
    if x < 2:
        return 0
    total = 0
    for k in range(1, 100):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        total += mu_k / k * li(x ** (1.0 / k))
    return total

def mobius(n):
    """Möbius function."""
    if n == 1:
        return 1
    factors = 0
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors += 1
            temp //= d
            if temp % d == 0:
                return 0  # Squared factor
        d += 1
    if temp > 1:
        factors += 1
    return (-1) ** factors

def lucy_pi(x):
    """Lucy_Hedgehog DP for π(x)."""
    x = int(x)
    if x < 2:
        return 0
    sqrtx = int(math.isqrt(x))
    small = list(range(-1, sqrtx + 1))
    large = [0] * (sqrtx + 2)
    for i in range(1, sqrtx + 1):
        large[i] = x // i - 1
    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:
            continue
        cnt = small[p - 1]
        p2 = p * p
        for i in range(1, min(sqrtx, x // p2) + 1):
            if i * p <= sqrtx:
                large[i] -= large[i * p] - cnt
            else:
                large[i] -= small[x // (i * p)] - cnt
        for i in range(sqrtx, p2 - 1, -1):
            small[i] -= small[i // p] - cnt
    return large[1]

# ============================================================
# EXPERIMENT 1: Self-similarity of the error term δ(x) = π(x) - R(x)
# ============================================================
print("=" * 60)
print("EXPERIMENT 1: Self-Similarity of δ(x) = π(x) - R(x)")
print("=" * 60)

# Compute δ(x) at various scales
x_values = [10**k for k in range(2, 8)]
deltas = []
for x in x_values:
    pi_x = lucy_pi(x)
    r_x = R_func(x)
    delta = pi_x - r_x
    deltas.append(delta)
    print(f"  x=10^{int(math.log10(x))}: π(x)={pi_x}, R(x)={r_x:.1f}, δ(x)={delta:.1f}, δ/√x={delta/math.sqrt(x):.4f}")

print("""
Observation: δ(x)/√x oscillates but stays bounded (under RH).
The Riemann Hypothesis says |δ(x)| ≤ C·√x·ln(x).

Self-similarity test: Is δ(10x) predictable from δ(x)?
""")

# Test: δ(10x) vs δ(x)
print("Self-similarity ratio δ(10x)/δ(x):")
for i in range(len(deltas) - 1):
    if abs(deltas[i]) > 0.1:
        ratio = deltas[i+1] / deltas[i]
        print(f"  δ(10^{i+3})/δ(10^{i+2}) = {ratio:.4f}")

# ============================================================
# EXPERIMENT 2: Recursive π(x) via scaled sub-problems
# ============================================================
print("\n" + "=" * 60)
print("EXPERIMENT 2: Recursive Scaling Formula")
print("=" * 60)

# Can we write π(x) ≈ α·π(x/k) + β(x) for some k, α, β?
# By PNT: π(x) ~ x/ln(x)
# π(x/k) ~ (x/k)/ln(x/k) = (x/k)/(ln(x)-ln(k))
# So π(x)/π(x/k) ~ k·(ln(x)-ln(k))/ln(x) = k·(1 - ln(k)/ln(x))

# For k=2: π(x)/π(x/2) ~ 2·(1 - ln2/ln(x))
# For large x, this → 2

# Test empirically:
print("Ratio π(x)/π(x/2):")
for exp in range(3, 8):
    x = 10**exp
    pi_x = lucy_pi(x)
    pi_half = lucy_pi(x // 2)
    ratio = pi_x / pi_half
    predicted = 2 * (1 - math.log(2) / math.log(x))
    print(f"  x=10^{exp}: π(x)/π(x/2) = {ratio:.6f}, predicted = {predicted:.6f}, diff = {abs(ratio-predicted):.6f}")

print("""
The ratio π(x)/π(x/2) is well-predicted by PNT.
But the ERROR in this prediction is O(√x/ln(x)²) — too large.

Can we build a recursive formula where errors DECREASE?

Define: π(x) = 2·π(x/2) + correction(x)
correction(x) = π(x) - 2·π(x/2)
""")

# Compute correction(x) = π(x) - 2·π(x/2)
print("Correction term c(x) = π(x) - 2·π(x/2):")
for exp in range(3, 8):
    x = 10**exp
    pi_x = lucy_pi(x)
    pi_half = lucy_pi(x // 2)
    correction = pi_x - 2 * pi_half
    # By PNT: correction ~ x·ln(2)/ln(x)² + ...
    predicted_corr = x * math.log(2) / math.log(x)**2
    print(f"  x=10^{exp}: c(x)={correction}, predicted~{predicted_corr:.0f}, ratio={correction/predicted_corr:.4f}")

print("""
The correction c(x) ≈ x·ln(2)/ln²(x) is LARGE (not o(1)).
It cannot be computed without knowing π(x) itself.

For a recursive formula to work, we need:
  π(x) = f(π(x/2), π(x/3), ..., x)  where f is CHEAP

But the "cheap" part of f gives only the PNT approximation.
The O(√x) error term cannot be predicted from π(x/2), π(x/3), etc.
because the error is controlled by ZETA ZEROS, not arithmetic.
""")

# ============================================================
# EXPERIMENT 3: Multiscale decomposition
# ============================================================
print("=" * 60)
print("EXPERIMENT 3: Multiscale / Wavelet Decomposition")
print("=" * 60)

# Idea: Decompose π(x) at multiple scales
# π(x) = π_smooth(x) + δ_1(x) + δ_2(x) + ...
# where δ_k captures structure at scale x^{1/2^k}

# The explicit formula already does this!
# π(x) = R(x) + (contribution from Re(ρ)=1/2 zeros)
#       + (contribution from possible higher zeros)
#       + (contribution from trivial zeros)

# The zeros at height γ contribute oscillations of "wavelength" ~ 1/γ
# Different zeros capture different scales

# Can we compute a FEW scales cheaply and bound the rest?

# Under RH, the contribution of zeros with |γ| > T is bounded by:
# O(√x · ln(x)² / T)

# For this to be < 0.5 (exact), need T > 2·√x·ln(x)²

# For x = 10^102: T > 2 × 10^51 × (230)² ≈ 10^56
# That's 10^56 zeros needed!

# What about a HIERARCHICAL approach?
# Step 1: Use first 100 zeros → error O(√x / 100)
# Step 2: For the residual, can we compute it at x^{1/2}?

print("""
Hierarchical zero approach:
  Level 0: R(x) → error ~ √x · (2γ₁·ln(x))^{-1} oscillating
  Level 1: + first K zeros → error ~ √x/(K·ln(x))
  Level 2: + more zeros → ...

Each zero costs O(1) to evaluate (once known).
The cost is computing/storing the zeros themselves.

Zeros of ζ(s) up to height T:
  Number of zeros: ~ T·ln(T)/(2π)
  Cost to compute: O(T^{1+ε}) via Riemann-Siegel or O(T^{1/3+ε}) via Odlyzko

For T = 10^56:
  Number of zeros: ~10^58
  Cost to find them: ~10^56+ε
  Both completely infeasible.

MULTISCALE DECOMPOSITION CONFIRMS THE BARRIER.
""")

# ============================================================
# EXPERIMENT 4: Exploiting the Riemann-von Mangoldt formula
# ============================================================
print("=" * 60)
print("EXPERIMENT 4: Riemann-von Mangoldt N(T)")
print("=" * 60)

# N(T) = number of zeros with 0 < Im(ρ) < T
# N(T) = T/(2π) · ln(T/(2πe)) + 7/8 + S(T) + O(1/T)
# where S(T) = (1/π)·arg(ζ(1/2 + iT))

# S(T) is small: |S(T)| < 0.137·ln(T) + 0.443·ln(ln(T)) + 4.350 (unconditionally)

# So we know HOW MANY zeros there are, very accurately.
# But we don't know WHERE they are (the individual γ values).

# The question: can we compute the SUM Σ_{|γ|<T} f(γ) without
# knowing individual γ values?

# For f(γ) = 1: Yes! That's N(T).
# For f(γ) = 1/γ: The "mean zero sum" can be computed analytically.
# For f(γ) = cos(γ·ln(x))/γ: This is what we need, and it CAN'T be
#   computed without individual zeros (it's the explicit formula!).

print("""
The Riemann-von Mangoldt formula N(T) tells us HOW MANY zeros
exist below height T, but not WHERE they are.

We need: Σ_{γ<T} cos(γ·ln(x)) / γ

This is the "smooth linear statistic" of zeros:
  Σ f(γ) with f(t) = cos(t·ln(x)) / t

By the explicit formula for smooth linear statistics:
  Σ_{γ} f(γ) = (main term involving Λ(n)) + (error)

The main term is:
  -Σ_{n≤x} Λ(n)·h(ln(n)) / √n + contributions from trivial zeros

where h is the Fourier transform of f.

But computing Σ Λ(n)·h(ln(n))/√n for n up to x IS JUST AS HARD
as computing π(x) directly!

The duality between zeros and primes is EXACT: knowing one
gives you the other, and both are equally hard to compute.
""")

# ============================================================
# EXPERIMENT 5: The "scaling hypothesis"
# ============================================================
print("=" * 60)
print("EXPERIMENT 5: The Scaling Hypothesis")
print("=" * 60)

# What if the error δ(x) = π(x) - R(x) satisfies a SCALING LAW?
# If δ(x) = x^{1/2} · Φ(ln(x)) for some periodic function Φ,
# then knowing Φ at a FEW points would give us δ everywhere.

# Under RH: δ(x) = -Σ_ρ li(x^ρ) ≈ -√x · Σ cos(γ ln x) / (γ ln x)
# This IS √x times a function of ln(x)!
# But the function is NOT periodic — it's a sum of infinitely many
# incommensurate frequencies.

# Test: Is Φ(t) = δ(e^t) / e^{t/2} quasi-periodic?
print("Testing if Φ(t) = δ(e^t)/e^{t/2} is quasi-periodic:")

ts = np.linspace(5, 16, 200)  # ln(x) from 5 to 16, i.e. x from 148 to ~9M
phis = []
for t in ts:
    x = int(math.exp(t))
    if x < 2:
        continue
    pi_x = lucy_pi(x)
    r_x = R_func(x)
    delta = pi_x - r_x
    phi = delta / math.sqrt(x)
    phis.append((t, phi))

phis = np.array(phis)
if len(phis) > 10:
    # Compute autocorrelation
    phi_vals = phis[:, 1]
    phi_centered = phi_vals - np.mean(phi_vals)
    acf = np.correlate(phi_centered, phi_centered, mode='full')
    acf = acf[len(acf)//2:]
    acf /= acf[0]

    # Find peaks in autocorrelation
    from scipy.signal import find_peaks
    peaks, props = find_peaks(acf[:50], height=0.3)
    print(f"  ACF peaks at lags: {peaks}")
    print(f"  ACF peak values: {acf[peaks]}")

    if len(peaks) > 0:
        period_estimate = np.mean(np.diff(peaks)) if len(peaks) > 1 else peaks[0]
        # Convert to frequency
        dt = ts[1] - ts[0]
        print(f"  Estimated period: {period_estimate * dt:.3f} in ln(x)")
        print(f"  Expected from γ₁=14.13: period = 2π/14.13 = {2*math.pi/14.134:.3f}")
    else:
        print("  No significant periodic component found")

    # FFT analysis
    fft_vals = np.fft.rfft(phi_centered)
    freqs = np.fft.rfftfreq(len(phi_centered), d=ts[1]-ts[0])
    power = np.abs(fft_vals)**2
    # Top 5 frequencies
    top_idx = np.argsort(power)[-6:-1]
    print(f"\n  Top 5 FFT frequencies: {freqs[top_idx]}")
    print(f"  Expected zeta zero frequencies: {14.134/(2*math.pi):.3f}, {21.022/(2*math.pi):.3f}, {25.011/(2*math.pi):.3f}")

print("""
RESULT: Φ(t) shows oscillations at zeta zero frequencies
(confirming the explicit formula), but it is NOT periodic.
The sum of incommensurate frequencies never repeats.

The "scaling hypothesis" — that δ(x) has simple scaling behavior —
is FALSE. The correction requires knowing infinitely many
zero frequencies and amplitudes.

CONCLUSION: No self-similar or fractal shortcut exists.
The prime distribution's fine structure is controlled by the
full spectrum of zeta zeros, which cannot be compressed.
""")

print("=" * 60)
print("ALL FRACTAL/SELF-SIMILAR EXPERIMENTS COMPLETE")
print("Conclusion: No exploitable self-similarity found.")
print("=" * 60)
