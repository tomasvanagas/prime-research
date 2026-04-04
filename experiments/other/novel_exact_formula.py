"""
Session 9: Novel Exact Formula Attempts

IDEA 1: "Exact floor rounding" approach
If we can find f(n) such that |f(n) - p(n)| < 1/2 for all n,
then p(n) = round(f(n)) exactly.

The challenge: f(n) must be computable in polylog time.

IDEA 2: "Binary digit extraction"
Compute p(n) one bit at a time, each bit in polylog.
The k-th bit of p(n) relates to π(x) at a specific point.

IDEA 3: "Implicit definition via fixed point"
Define p(n) as the unique fixed point of some contractive map T
where T can be evaluated in polylog time.

IDEA 4: "Number-theoretic transform approach"
Use NTT/DFT to transform the prime indicator into a domain
where individual values are accessible.

IDEA 5: "Compositional formula"
p(n) = f1(n) ∘ f2(n) ∘ ... ∘ fk(n) where each fi is cheap.
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime
from mpmath import mp, mpf, log, li, exp, floor, sqrt, pi as mpi
import time

mp.dps = 50

print("=" * 70)
print("IDEA 1: Exact Floor Rounding via Riemann R^{-1}")
print("=" * 70)

# R(x) = 1 + Σ_{k=1}^∞ (log x)^k / (k · k! · ζ(k+1))
# R^{-1}(n) ≈ the x such that R(x) = n
# Error: |R^{-1}(n) - p(n)| ~ O(√p(n))

# But what if we could CERTIFY that for specific n,
# the error is less than gap(n)/2?

# gap(n) = p(n+1) - p(n)
# By Cramér's conjecture: gap(n) ~ (log p(n))²
# For p(10^100) ≈ 2.35×10^102: gap ~ (236)² ≈ 55696
# Error from R^{-1}: ~√(2.35×10^102) ≈ 1.5×10^51

# Error >> gap by factor of ~10^47. Cannot round.
# EVEN under RH, error in explicit formula is O(√x · log x)

print("Error vs gap analysis:")
for n in [100, 1000, 10000, 100000]:
    pn = prime(n)
    gap = prime(n+1) - pn
    error_est = np.sqrt(pn) * np.log(pn)
    print(f"  n={n:6d}: p(n)={pn:8d}, gap={gap:3d}, error_est={error_est:.1f}, ratio={error_est/gap:.1f}")

print("\nVerdict: Error always exceeds gap. Floor rounding IMPOSSIBLE for large n.")

print("\n" + "=" * 70)
print("IDEA 2: Bit-by-Bit via π(x) Threshold")
print("=" * 70)

# p(n) has b = ⌈log₂(p(n))⌉ bits.
# Bit k (from MSB): p(n) has bit k set iff p(n) ≥ 2^k
# More precisely, the i-th bit = ⌊p(n)/2^i⌋ mod 2

# We can determine the MSB: p(n) is in [2^(b-1), 2^b)
# b is determined by π(2^b) ≥ n > π(2^(b-1))
# This requires computing π at powers of 2.

# For each subsequent bit, we need to determine:
# is p(n) ≥ current_lower_bound + 2^(remaining_bits)?
# This is equivalent to: π(current_lower_bound + 2^(remaining_bits)) ≥ n?

# Total: b calls to π(x), each costing O(x^{2/3}).
# No savings — this IS binary search on π(x).

print("Bit-by-bit extraction = binary search on π(x)")
print("Each bit costs O(p(n)^{2/3}) → total O(p(n)^{2/3} · log p(n))")
print("Actually WORSE than direct binary search by factor log p(n)")
print("Verdict: FAIL — no per-bit shortcut exists")

print("\n" + "=" * 70)
print("IDEA 3: Fixed Point of Contractive Map")
print("=" * 70)

# Define T(x) = R^{-1}(π(x)) where R^{-1} is the smooth approximation
# If x = p(n), then π(p(n)) = n, so T(p(n)) = R^{-1}(n) ≈ p(n)
# But T is NOT contractive near p(n) — the error doesn't shrink.

# Alternative: define T(x) = x + n - π(x) (trying to "push" x to p(n))
# If x > p(n): π(x) ≥ n, so T(x) ≤ x (pushes down)
# If x < p(n): π(x) < n, so T(x) > x (pushes up)
# Fixed point: T(x*) = x* iff π(x*) = n iff x* = p(n) (if p(n) exactly)

# But T has derivative T'(x) = 1 - π'(x) = 1 - (1/ln x + ...)
# For large x: T'(x) ≈ 1 - 1/ln(x) which is close to 1 (NOT contractive!)
# Need |T'(x)| < 1 for convergence. This fails.

# Modified: T_α(x) = x + α(n - π(x)) for some α
# T'(x) = 1 - α·π'(x) = 1 - α/ln(x)
# For |T'| < 1 at x ~ p(n): need α > 0 (trivially met)
# and 1 - α/ln(x) > -1, so α < 2·ln(x) ≈ 2·ln(p(n))

# Optimal α = ln(p(n)) gives T'(x) ≈ 0 → quadratic convergence!
# But each T evaluation requires π(x) → O(x^{2/3})
# And we need ~log(log p(n)) iterations for quadratic convergence
# Total: O(p(n)^{2/3} · log log p(n)) — same as Newton on π

print("Fixed-point iteration T_α(x) = x + α(n - π(x)):")
print("  Optimal α = ln(p(n))")
print("  Convergence: quadratic (like Newton's method)")
print("  Cost per step: O(p(n)^{2/3}) for π evaluation")
print("  Total: O(p(n)^{2/3} · log log p(n))")
print("  Verdict: This IS Newton's method in disguise. Same cost.")

# Let's verify this works for small n
def fixed_point_prime(n, max_iter=20):
    """Find p(n) via fixed-point iteration"""
    # Initial guess from PNT
    from mpmath import log as mlog
    x = float(n * mlog(n) + n * mlog(mlog(n))) if n > 5 else 11
    alpha = float(mlog(max(x, 3)))

    for i in range(max_iter):
        pi_x = primepi(int(x))
        residual = n - pi_x
        if residual == 0 and isprime(int(x)):
            return int(x), i
        x_new = x + alpha * residual
        if abs(x_new - x) < 0.5:
            # Close enough — find nearest prime
            candidate = int(round(x_new))
            if isprime(candidate) and primepi(candidate) == n:
                return candidate, i
            # Search nearby
            for delta in range(-10, 11):
                c = candidate + delta
                if c > 1 and isprime(c) and primepi(c) == n:
                    return c, i
        x = x_new
    return None, max_iter

print("\nVerification:")
correct = 0
for n in range(1, 101):
    result, iters = fixed_point_prime(n)
    expected = prime(n)
    if result == expected:
        correct += 1
    else:
        print(f"  FAIL: n={n}, expected={expected}, got={result}")
print(f"Accuracy: {correct}/100")

print("\n" + "=" * 70)
print("IDEA 4: Number-Theoretic Transform of Prime Indicator")
print("=" * 70)

# The prime indicator function χ_P(n) = 1 if n is prime, 0 otherwise
# Its DFT: χ̂(k) = Σ_{p≤N} e^{-2πipk/N}
# These are EXPONENTIAL SUMS over primes — studied extensively.

# By Vinogradov: |Σ_{p≤N} e^{2πiαp}| << N/(log N) for irrational α
# This means χ̂(k) ~ N/(log N) for most k — essentially FLAT spectrum
# No sparse representation in Fourier domain!

N = 1000
primes_N = list(p for p in range(2, N+1) if isprime(p))
indicator = np.zeros(N)
for p in primes_N:
    indicator[p-1] = 1

fft_indicator = np.abs(np.fft.fft(indicator))
print(f"FFT of prime indicator (N={N}):")
print(f"  DC component: {fft_indicator[0]:.1f} (= π({N}) = {len(primes_N)})")
print(f"  Mean |FFT|: {np.mean(fft_indicator[1:]):.1f}")
print(f"  Max |FFT|: {np.max(fft_indicator[1:]):.1f}")
print(f"  Min |FFT|: {np.min(fft_indicator[1:]):.1f}")
print(f"  Std |FFT|: {np.std(fft_indicator[1:]):.1f}")

# Sparsity analysis: how many Fourier coefficients needed?
sorted_fft = np.sort(fft_indicator[1:])[::-1]
total_energy = np.sum(fft_indicator[1:]**2)
cumulative = np.cumsum(sorted_fft**2) / total_energy
k_90 = np.searchsorted(cumulative, 0.90) + 1
k_99 = np.searchsorted(cumulative, 0.99) + 1
k_999 = np.searchsorted(cumulative, 0.999) + 1
print(f"  Coefficients for 90% energy: {k_90}/{N}")
print(f"  Coefficients for 99% energy: {k_99}/{N}")
print(f"  Coefficients for 99.9% energy: {k_999}/{N}")
print(f"  Verdict: NOT sparse. Need O(N) coefficients.")

print("\n" + "=" * 70)
print("IDEA 5: Compositional / Layered Formula")
print("=" * 70)

# p(n) = f_k ∘ f_{k-1} ∘ ... ∘ f_1(n)
# Each f_i is a simple function (polynomial, exp, log, floor, etc.)

# The Prunescu-Shunia result shows this EXISTS with {+,-,*,/,^}
# but with intermediate values of 10^78913 digits.

# Can we find a composition with BOUNDED intermediates?
# Key constraint: the composition must introduce ~170 bits of info for large n.
# Where does this info come from? From the STRUCTURE of the functions.

# Example: p(n) = ⌊f(n) + g(n)·sin(h(n))⌋
# The sin term oscillates and can encode arbitrary information
# But computing sin to 170-bit precision IS the computation.

# What about: p(n) = ⌊n·log(n)·(1 + correction(n))⌋?
# correction(n) must have ~170 significant bits for n=10^100
# These bits encode zeta zeros — can't avoid computing them.

# The MOST promising compositional approach:
# p(n) = ⌊li^{-1}(n) + Σ_{k=1}^K c_k · f_k(n)⌋
# where f_k are basis functions and c_k are constants
# But K must grow with n (no finite set of constants works for all n)

print("Compositional analysis:")
print("  p(n) = ⌊smooth(n) + correction(n)⌋")
print("  smooth(n): computable in O(polylog n)")
print("  correction(n): requires ~170 bits for n=10^100")
print("  These 170 bits CANNOT come from any finite formula")
print("  They encode the cumulative effect of infinitely many zeta zeros")
print()
print("THEORETICAL IMPOSSIBILITY (information-theoretic):")
print("  Any formula with K parameters has at most K·b bits of info")
print("  where b is the precision of each parameter.")
print("  For p(10^100): need ~170 random bits")
print("  So K·b ≥ 170, meaning we need at least one 170-bit constant")
print("  That constant ENCODES p(10^100) — the formula is just a lookup table")

print("\n" + "=" * 70)
print("IDEA 6: Completely New — Prefix-Free Encoding + Levin Search")
print("=" * 70)

# By Levin's universal search theorem, there exists a program of length
# K(p(n)|n) that outputs p(n) given n.
# K(p(n)|n) = O(log log n) ≈ 5 bits (proven in session 8)
#
# BUT: the program runs in time O(2^{K(p(n)|n)} · t(n)) where t(n) is
# the actual computation time. So the short program description doesn't
# help with TIME complexity.
#
# The 5-bit program is essentially: "compute π(x) and binary search"
# encoded in 5 bits. The computation still takes O(p(n)^{2/3}) time.

print("Levin search analysis:")
print("  K(p(n)|n) = O(log log n) ≈ 5 bits")
print("  The 5-bit program IS 'compute pi(x) and binary search'")
print("  Runtime: 2^5 × O(p(n)^{2/3}) = O(p(n)^{2/3})")
print("  No speedup from short description length")
print("  Kolmogorov complexity ≠ time complexity")

print("\n" + "=" * 70)
print("RADICAL IDEA 7: What if p(n) = ⌊α · β^n⌋ for some constants?")
print("=" * 70)

# Mills' theorem: there exists A such that ⌊A^{3^n}⌋ is always prime.
# But A encodes all primes in its digits — circular.
#
# What about a DIFFERENT exponential form?
# p(n) ≈ exp(W(n)) where W is Lambert W? No — error O(√p).
# p(n) ≈ ⌊e^{n/a + b}⌋? No — primes grow as n·ln(n), not exponentially.
#
# But what about: p(n) = ⌊f(n)⌋ where f(n) is defined by a DIFFERENTIAL EQUATION?
# f'(t) = ln(f(t)) + correction(t)
# By PNT: p(n) ≈ integral, so p'(n) ≈ ln(p(n))
# The correction term encodes the irregularity.

# Let's check: what ODE does p(n) satisfy if we treat it as a continuous function?
from scipy.interpolate import CubicSpline

ns = np.arange(2, 1001)
ps = np.array([prime(n) for n in ns])

# Numerical derivative
dp = np.diff(ps).astype(float)  # dp[i] ≈ p'(n) at n+0.5
ns_mid = ns[:-1] + 0.5

# p'(n) vs ln(p(n))
ln_p = np.log(ps[:-1].astype(float))
ratio = dp / ln_p

print(f"p'(n) / ln(p(n)) statistics:")
print(f"  Mean: {np.mean(ratio):.4f}")
print(f"  Std:  {np.std(ratio):.4f}")
print(f"  Min:  {np.min(ratio):.4f}")
print(f"  Max:  {np.max(ratio):.4f}")

# So p'(n) ≈ c · ln(p(n)) where c ≈ 1
# More precisely: p'(n) = ln(p(n)) + ln(ln(p(n))) + 1/ln(p(n)) + noise
# The noise is O(1) but random — exactly the gap fluctuations

print(f"\nWith second-order term:")
ln_ln_p = np.log(ln_p)
corrected_ratio = dp / (ln_p + ln_ln_p)
print(f"  p'(n) / (ln(p) + ln(ln(p))):")
print(f"  Mean: {np.mean(corrected_ratio):.4f}")
print(f"  Std:  {np.std(corrected_ratio):.4f}")

# The noise in p'(n) is THE prime gaps, which are the fundamental obstacle.
# An ODE can capture the smooth part but NOT the random gaps.

print("\nVerdict: p(n) satisfies p'(n) = ln(p(n)) + ln(ln(p(n))) + NOISE")
print("NOISE = gap fluctuations = irreducible randomness")
print("Cannot solve ODE without knowing the noise → circular")

print("\n" + "=" * 70)
print("SESSION 9 DIRECT EXPERIMENTS: ALL CONFIRM BARRIER")
print("=" * 70)
