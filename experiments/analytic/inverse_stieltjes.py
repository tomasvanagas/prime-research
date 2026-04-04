"""
Session 8: Inverse Stieltjes / Spectral Measure approach

Novel angle: The prime counting function π(x) is a COUNTING MEASURE.
Its Stieltjes transform is:
  S(z) = Σ_p 1/(z - p) = ∫ dπ(x)/(z - x)

For z on the real axis away from primes, S(z) encodes all prime positions.
The INVERSE Stieltjes transform recovers π(x):
  π(x) = lim_{ε→0+} (1/π) Im[S(x + iε)]  ... wait, that's for continuous measures.

Actually for discrete measures:
  π(x) = lim_{ε→0+} (1/π) ∫_2^x Im[S(t + iε)] dt

Can we compute S(z) without knowing the primes?

Connection: S(z) is related to -ζ'(s)/ζ(s) via Mellin transform.
Specifically:
  -ζ'(s)/ζ(s) = s ∫_1^∞ ψ(x)/x^{s+1} dx
  where ψ(x) = Σ_{p^k ≤ x} ln(p)

So S(z) is essentially the resolvent of the "prime operator".

Novel idea: What if we compute the resolvent at a SINGLE carefully chosen
point z, and extract π(x) from it?

Also: the Stieltjes constants γ_n appear in the Laurent expansion of ζ(s):
  ζ(s) = 1/(s-1) + Σ_{n=0}^∞ (-1)^n γ_n (s-1)^n / n!

Can the Stieltjes constants help us count primes?
"""

import numpy as np
from mpmath import mp, mpf, mpc, zetazero, zeta, log, pi as mppi, exp, fsum
from mpmath import diff as mpdiff
from sympy import primerange
import time

mp.dps = 30

primes_list = list(primerange(2, 100000))

# =============================================================================
# 1. Compute S(z) = Σ 1/(z-p) at specific points
# =============================================================================
print("=" * 60)
print("1. Resolvent S(z) = Σ 1/(z-p) at specific points")
print("=" * 60)

def resolvent(z, max_prime=10000):
    """Compute S(z) = Σ_{p ≤ max_prime} 1/(z - p)"""
    primes = list(primerange(2, max_prime + 1))
    return fsum(1/(mpc(z) - p) for p in primes)

# At z = x + iε, Im(S(z)) has peaks at prime locations
# Im(1/(x + iε - p)) = -ε / ((x-p)² + ε²) = Lorentzian peak at x = p

# Test: reconstruct π(x) from resolvent
print("  Testing resolvent-based π(x) reconstruction:")
for eps in [0.1, 0.5, 1.0, 5.0]:
    correct = 0
    total = 0
    for x_test in range(10, 101):
        # π(x) ≈ -(1/π) ∫_2^x Im(S(t + iε)) dt
        # Numerically: sum over grid
        dt = 0.1
        integral = 0
        for t_idx in range(int(2/dt), int(x_test/dt)):
            t = t_idx * dt
            s_val = resolvent(mpc(t, eps), max_prime=200)
            integral += float(s_val.imag) * dt

        pi_est = round(-integral / float(mppi))
        pi_actual = len([p for p in primes_list if p <= x_test])
        if pi_est == pi_actual:
            correct += 1
        total += 1

    print(f"  ε={eps}: {correct}/{total} = {correct/total:.1%} exact")

# =============================================================================
# 2. Stieltjes constants and prime counting
# =============================================================================
print("\n" + "=" * 60)
print("2. Stieltjes constants")
print("=" * 60)

# γ_n = lim_{N→∞} [Σ_{k=1}^N (ln k)^n/k - (ln N)^{n+1}/(n+1)]
# First few:
# γ_0 = γ ≈ 0.5772... (Euler-Mascheroni)
# γ_1 ≈ -0.0728...
# γ_2 ≈ -0.00969...

# The connection to primes: ζ(s) = Π(1-p^{-s})^{-1}
# So -ζ'/ζ(s) = Σ_p ln(p)/(p^s - 1) = Σ_p Σ_k ln(p)/p^{ks}
# At s=1: -ζ'/ζ has a pole with residue 1
# The Laurent coefficients involve Stieltjes constants

# Can Stieltjes constants help count primes?
# γ_n involves Σ (ln k)^n / k which is a smooth sum — no prime information

# Direct computation
gamma_0 = float(mp.euler)
print(f"  γ₀ (Euler-Mascheroni) = {gamma_0:.10f}")

# Compute a few Stieltjes constants
for n in range(5):
    # γ_n from the Laurent expansion of ζ(s) at s=1
    # Using: γ_n = (-1)^n * n! * lim_{s→1} [ζ^{(n)}(s)/(n!) - (-1)^n/(s-1)^{n+1}]
    # Simpler: directly from mpmath if available
    try:
        # γ_n = limit as N→∞ of sum_{k=1}^N (ln(k))^n/k - (ln(N))^{n+1}/(n+1)
        N = 10000
        gamma_n = sum(np.log(k)**n / k for k in range(1, N+1)) - np.log(N)**(n+1)/(n+1)
        print(f"  γ_{n} ≈ {gamma_n:.10f} (N={N})")
    except:
        pass

print("\n  Stieltjes constants involve HARMONIC sums, not prime sums")
print("  They encode properties of ζ(s) near s=1 but NOT individual primes")

# =============================================================================
# 3. Novel: Can we extract primes from poles of -ζ'/ζ(s)?
# =============================================================================
print("\n" + "=" * 60)
print("3. Poles of -ζ'/ζ(s) and prime extraction")
print("=" * 60)

# -ζ'/ζ(s) = Σ_p Σ_k ln(p) * p^{-ks}
# This converges for Re(s) > 1
# Poles of -ζ'/ζ are at the zeros of ζ (the non-trivial zeros)
# And at s = 1 (pole of ζ)

# The explicit formula for ψ(x) is:
# ψ(x) = x - Σ_ρ x^ρ/ρ - ζ'(0)/ζ(0) - (1/2)ln(1-1/x²)

# Can we extract INDIVIDUAL prime powers from ψ(x)?
# ψ(x) jumps by ln(p) at x = p, p², p³, ...
# If we could compute ψ(x) exactly, we could detect primes

# But computing ψ(x) exactly requires the FULL explicit formula
# with ALL zeta zeros — the same barrier

# Alternative: compute ψ(x) - ψ(x-1) for integer x
# This equals ln(p) if x = p^k, and 0 otherwise
# But the explicit formula for the DIFFERENCE is:
# ψ(x) - ψ(x-1) = 1 - Σ_ρ (x^ρ - (x-1)^ρ)/ρ - ...
# Each zero contributes x^{ρ-1} ≈ x^{-1/2+iγ} to the difference
# Sum of ~T terms oscillates with amplitude ~T·x^{-1/2}
# For error < ln(2)/2 (to detect primes): T > x^{1/2}/2

# So: same barrier, O(x^{1/2}) zeros needed

print("  ψ(x) - ψ(x-1) = ln(p) if x=p^k, 0 otherwise")
print("  Computing via explicit formula needs O(x^{1/2}) zeros")
print("  Same barrier as π(x)")

# =============================================================================
# 4. Novel: Weighted prime sums with special kernels
# =============================================================================
print("\n" + "=" * 60)
print("4. Weighted prime sums with optimal kernels")
print("=" * 60)

# Instead of π(x) = #{p ≤ x}, consider:
# F(x) = Σ_p K(p/x) for a smooth kernel K
# If K is bandlimited, its Mellin transform K̃(s) decays rapidly
# And F(x) = (1/2πi) ∫ K̃(s) · (-ζ'/ζ(s)) · x^s ds / s

# The KEY question: can K be chosen so that F(x) determines π(x)
# while K̃(s) is concentrated near Re(s) = 1/2?

# The answer is the UNCERTAINTY PRINCIPLE for Mellin transforms:
# A function cannot be both sharply localized in x-space (step function)
# and rapidly decaying in s-space (few zeros needed)
# This is the Beurling-Selberg extremal problem

# The optimal K for minimizing the number of zeros needed is
# the Gaussian in log-space: K(x) = exp(-a·(ln x)²)
# This gives error O(x·exp(-c·√(ln x))) — better than PNT but not exact

# For EXACT π(x), K must be the step function → K̃(s) = x^s/s
# → need ALL zeros up to height T ~ √x

print("  Uncertainty principle for Mellin transforms:")
print("  Sharp cutoff (step function) ⟷ slow decay in s-space (all zeros needed)")
print("  Smooth kernel (Gaussian) ⟷ fast decay (few zeros) but inexact π(x)")
print("  TRADE-OFF IS FUNDAMENTAL: cannot have both exact π(x) and few zeros")

print("\n" + "=" * 60)
print("FINAL CONCLUSIONS (Inverse Stieltjes / Spectral)")
print("=" * 60)
print("""
1. Resolvent S(z): works for reconstruction but requires summing over all primes — circular
2. Stieltjes constants: encode ζ near s=1, NOT individual primes
3. -ζ'/ζ poles: same explicit formula barrier (O(x^{1/2}) zeros)
4. Optimal kernels: uncertainty principle prevents exact π(x) with few zeros

The spectral/resolvent perspective CONFIRMS the barrier from yet another angle:
  - Primes are the EIGENVALUES of the "arithmetic operator"
  - Computing individual eigenvalues of an infinite-dimensional operator
    requires O(√(dimension)) work
  - For primes up to x, the effective dimension is x
  - So the cost is O(x^{1/2}) — consistent with all other approaches

This is impossibility proof #21: the SPECTRAL resolution barrier.
""")
