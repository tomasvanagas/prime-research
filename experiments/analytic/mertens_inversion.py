"""
Session 4: Mertens Function & Novel Möbius Inversion Approach

Key insight: π(x) and M(x) = Σ_{n≤x} μ(n) are deeply connected.
If we could compute M(x) in O(polylog), we could potentially extract π(x).

Also exploring: Can we use the INVERSE relationship?
  Σ_{d|n} μ(d) = [n=1]  (Möbius inversion)
  π(x) = Σ_{n=1}^{x} Σ_{d|n} Λ(d)/ln(d)  ... circular

New angle: The Meissel-Mertens constant and prime reciprocal sums
  Σ_{p≤x} 1/p = ln(ln(x)) + M₁ + O(1/ln(x))

  Can we use MULTIPLE such identities at different scales to triangulate π(x)?
"""

import math
import time
from functools import lru_cache

# First, let's explore: what if we could compute π(x) from
# MULTIPLE analytic functions evaluated at different points?

# The idea: We know many "sum over primes" identities:
# S_k(x) = Σ_{p≤x} p^k  for k = -1, 0, 1, 2, ...
# S_0(x) = π(x)
# S_{-1}(x) = Σ 1/p = ln(ln(x)) + M₁ + error
# S_1(x) = Σ p ≈ x²/(2 ln x) + ...
# S_k(x) has asymptotic formulas

# If we know S_k(x) exactly for k = -1, -2, ..., -K,
# can we recover S_0(x) = π(x)?

# The prime zeta function: P(s) = Σ_p p^{-s}
# P(s) = Σ_{k=1}^∞ μ(k)/k · ln ζ(ks)
# This converges for Re(s) > 1

# Key insight: P(s) for s > 1 is computable from ζ(s) which is computable!
# P(s) encodes ALL primes. Can we extract individual primes?

# Let's think about this differently using generating functions.
# Define f(z) = Σ_p z^p = z^2 + z^3 + z^5 + z^7 + z^11 + ...
# Then f(z) for |z| < 1 encodes all primes
# f(z) = Σ_{n=1}^∞ z^{p_n}
# The n-th prime is the position of the n-th non-zero coefficient

# But evaluating f(z) requires knowing primes...
# UNLESS we can express f(z) in terms of known functions.

# f(z) = Σ_{n=2}^∞ χ_prime(n) z^n
# where χ_prime is the prime indicator function

# Von Mangoldt: χ_prime(n) = (1/ln n) Σ_{d|n} μ(n/d) ln(d)... still circular

# NEW APPROACH: Perron's formula applied to prime zeta function
# P(s) = Σ_p p^{-s}
# By Perron's formula:
# π(x) = (1/2πi) ∫_{c-i∞}^{c+i∞} P(s) x^s / s ds
#
# Wait... P(s) = Σ_p p^{-s}, so
# (1/2πi) ∫ P(s) x^s/s ds = Σ_p (1/2πi) ∫ (x/p)^s / s ds = Σ_p [p ≤ x] = π(x)
#
# But P(s) has a natural boundary at Re(s) = 0, so the contour can't be shifted!
#
# However: P(s) = ln ζ(s) - Σ_{k=2}^∞ P(ks)/k
#         P(s) = ln ζ(s) - P(2s)/2 - P(3s)/3 - ...
#
# For Re(s) > 1, the series P(ks) for k≥2 converge VERY fast.
# P(2s) for Re(s) > 1 means evaluating at Re > 2, where the series has < 0.35 terms effectively
# So P(s) ≈ ln ζ(s) - P(2s)/2 with exponentially small corrections

# Let me try a completely different approach:
# INCLUSION-EXCLUSION ON SMOOTH NUMBERS

# Let S(x,y) = #{n ≤ x : all prime factors of n > y}
# Then π(x) - π(y) + 1 = S(x,y) for y = √x (Legendre's identity!)
# More precisely: π(x) = S(x,√x) + π(√x) - 1

# S(x,y) satisfies Buchstab's identity:
# S(x,y) = S(x,z) - Σ_{y≤p<z} S(x/p, p)  for z > y

# This is recursive... but what if we can compute S(x,y) analytically?
# The Dickman function ρ(u) gives: S(x, x^{1/u}) ~ x ρ(u) / ln(x^{1/u})
# ρ(u) satisfies: u ρ'(u) = -ρ(u-1) for u > 1, ρ(u) = 1 for 0 ≤ u ≤ 1

# The Dickman function is COMPUTABLE in O(polylog) time!
# But S(x,y) ~ x ρ(u) is only asymptotic...

# Let me actually TEST: how good is the Dickman approximation?

def dickman_rho(u, terms=50):
    """Compute Dickman's function ρ(u) using series expansion."""
    if u <= 1:
        return 1.0
    if u <= 2:
        return 1.0 - math.log(u)
    # For u > 2, use numerical ODE or recursion
    # Simple: use piecewise for small u
    if u <= 3:
        # ρ(u) = 1 - ln(u) + ∫₂ᵘ ln(t-1)/t dt for 2 < u ≤ 3
        from scipy.integrate import quad
        val, _ = quad(lambda t: math.log(t-1)/t, 2, u)
        return 1 - math.log(u) + val
    # For larger u, ρ(u) ≈ exp(-u(ln u + ln ln u - 1))
    return math.exp(-u * (math.log(u) + math.log(math.log(u)) - 1))


# Actually let me focus on something more promising.
# APPROACH: Use the EXPLICIT formula with a FINITE number of zeros
# but with a SMART correction term.

# The explicit formula: π(x) = li(x) - Σ_ρ li(x^ρ) - ln 2 + ∫_x^∞ dt/(t(t²-1)ln t)
#
# With K zeros: π_K(x) = li(x) - Σ_{|γ|≤T} li(x^ρ) - ln 2 + ...
# Error: |π(x) - π_K(x)| ≤ C x ln²x / T
#
# For x ~ 10^102, T = 10^48 needed. But...
#
# What if instead of trying to make the error < 1,
# we use the PARTIAL sum to get π(x) mod M for several moduli M,
# then use CRT to reconstruct π(x)?
#
# The error in π_K(x) is bounded but OSCILLATORY.
# If we evaluate at x and x+1, the ERROR changes slowly while π changes by 0 or 1.
# So: π(x+1) - π(x) = (π_K(x+1) - π_K(x)) + (small error change)
# The primality of x+1 is encoded in whether π changes!

# Actually this is still O(error) accuracy...

# Let me try: MODULAR PRIME COUNTING
# Compute π(x) mod p for small primes p, using properties of primes in APs

print("=" * 70)
print("SESSION 4: Novel Approaches to Prime Counting")
print("=" * 70)

# APPROACH 1: Can we compute π(x) mod m for small m?
#
# π(x) mod 2: parity of π(x)
# By Rubinstein-Sarnak (1994), the "parity" of π(x) relates to zeros of L-functions
# li(x) - π(x) changes sign infinitely often (Littlewood)
#
# But computing π(x) mod 2 exactly is AS HARD as computing π(x) exactly
# (because if you know π(x) mod 2 for all x, you know which x are prime)

# APPROACH 2: PRIME SUMS identity
# Σ_{p≤x} 1 = π(x)
# Σ_{p≤x} p = prime sum function
# Both are related via Abel summation:
# Σ_{p≤x} f(p) = f(x)π(x) - ∫₂ˣ f'(t)π(t) dt
#
# What if we use MULTIPLE such identities to set up a system of equations?
# For f(p) = p^k, k = 0, 1, 2, ..., K:
# S_k(x) = x^k π(x) - k ∫₂ˣ t^{k-1} π(t) dt
#
# Each S_k(x) can be approximated by the explicit formula...
# But they all require knowing π(t) for t ≤ x, which is circular.

# APPROACH 3: THE FIBONACCI-PRIME ANALOGY
# Fibonacci: F(n) = round(φ^n / √5) - EXACT from a simple formula
# Can we find an analogous formula for primes?
#
# Fibonacci satisfies a LINEAR recurrence with CONSTANT coefficients.
# Primes don't satisfy any linear recurrence (they grow too irregularly).
#
# But what about a NON-LINEAR recurrence?
# p(n+1) = p(n) + g(n), where g(n) is the gap
# g(n) ~ ln(p(n)) on average, but fluctuates wildly
#
# Cramér's model: g(n) ~ Poisson(ln p(n))
# This makes p(n) a sum of Poisson random variables... inherently random

# APPROACH 4: INFORMATION-THEORETIC ENCODING
# p(n) has about log₂(p(n)) ≈ n ln(n) ln(2) bits of information
# R^{-1}(n) captures about half these bits
# The remaining bits come from the oscillatory part of the explicit formula
#
# What if we could COMPRESS the oscillatory information?
# The zeta zeros are deterministic - they can be computed.
# So the "randomness" in δ(n) is actually deterministic!
# It's pseudo-random, not truly random.
#
# Can we find a PATTERN in the pseudo-randomness?

# Let me actually compute and analyze δ(n) for small n using known primes

def primes_up_to(n):
    """Simple sieve for reference primes."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = False
    return [i for i in range(2, n+1) if sieve[i]]

def li(x):
    """Logarithmic integral using series."""
    if x <= 1:
        return 0
    lnx = math.log(x)
    # Ramanujan's series for li(x)
    gamma = 0.5772156649015329
    result = gamma + math.log(abs(lnx))
    term = lnx
    for k in range(1, 100):
        result += term / (k * math.factorial(k))
        term *= lnx
        if abs(term / (k * math.factorial(k))) < 1e-15:
            break
    return result

def R_inverse(n):
    """Riemann R inverse: solve R(x) = n via Newton's method."""
    if n <= 1:
        return 2
    # Initial guess from PNT
    x = n * math.log(n)
    if n > 5:
        x = n * (math.log(n) + math.log(math.log(n)))

    for _ in range(100):
        if x <= 1:
            x = 2.0
        # R(x) = Σ μ(k)/k * li(x^{1/k})
        rx = 0
        for k in range(1, 50):
            mu_k = mobius(k)
            if mu_k == 0:
                continue
            xk = x ** (1.0/k)
            if xk < 1.001:
                break
            rx += mu_k / k * li(xk)

        # R'(x) = 1/(x ln x) approximately (dominant term)
        if x <= 1:
            x = 2.0
        rprime = 1.0 / (x * math.log(x))

        dx = (n - rx) / rprime
        x += dx
        if x < 2:
            x = 2.0
        if abs(dx) < 0.5:
            break

    return x

_mobius_cache = {}
def mobius(n):
    """Möbius function."""
    if n in _mobius_cache:
        return _mobius_cache[n]
    if n == 1:
        return 1
    # Factor n
    result = 1
    temp = n
    for p in range(2, int(n**0.5) + 1):
        if temp % p == 0:
            temp //= p
            if temp % p == 0:
                _mobius_cache[n] = 0
                return 0
            result *= -1
    if temp > 1:
        result *= -1
    _mobius_cache[n] = result
    return result

# Get reference primes
print("\nGenerating reference primes...")
primes = primes_up_to(200000)
print(f"Generated {len(primes)} primes up to {primes[-1]}")

# Compute δ(n) = p(n) - R^{-1}(n) for various n
print("\n--- Analyzing δ(n) = p(n) - R^{-1}(n) ---")
deltas = []
for n in range(2, min(len(primes), 2001)):
    ri = R_inverse(n)
    delta = primes[n-1] - ri
    deltas.append(delta)

print(f"Computed {len(deltas)} delta values")
print(f"Mean delta: {sum(deltas)/len(deltas):.4f}")
print(f"Std delta: {(sum(d**2 for d in deltas)/len(deltas) - (sum(deltas)/len(deltas))**2)**0.5:.4f}")
print(f"Max |delta|: {max(abs(d) for d in deltas):.4f}")

# APPROACH 5: DIGIT-BY-DIGIT COMPUTATION
# What if we can compute p(n) one digit at a time?
#
# Let p(n) = d_1 * 10^k + d_2 * 10^{k-1} + ... + d_{k+1}
# To find d_1: we need to know which interval [a*10^k, (a+1)*10^k) contains p(n)
# This requires: π((a+1)*10^k) - π(a*10^k) primes in the interval
# And summing to find n
#
# But computing π(x) is the bottleneck!
# UNLESS... we can compute π(x) mod 10 or something using the explicit formula

# APPROACH 6: ALGEBRAIC INDEPENDENCE OF ZETA ZEROS
# The first few zeta zeros ρ_k = 1/2 + iγ_k contribute:
# -li(x^ρ_k) ≈ -sqrt(x)/ln(x) * 2*cos(γ_k * ln(x))
#
# These are OSCILLATIONS with incommensurate frequencies γ_k
# By Kronecker's theorem, they're dense on the torus
# This means the correction is equidistributed... and unpredictable from n alone
#
# BUT: it's deterministic from x = p(n)!
# If we KNEW x, we could compute the correction.
# Chicken-and-egg: we need x to compute the correction, need the correction to find x.

# APPROACH 7: SELF-CONSISTENT ITERATION
# Start with x₀ = R^{-1}(n)
# Compute correction: c₁ = -Σ_{k=1}^{K} li(x₀^ρ_k)  [using first K zeros]
# Update: x₁ = R^{-1}(n + c₁)  ... no, this doesn't make sense
#
# Better: x₁ = x₀ + (n - π̃(x₀)) * ln(x₀)
# where π̃(x) = li(x) - Σ_{k=1}^K li(x^ρ_k) - ln2 + ...
# This is Newton's method on the explicit formula!
#
# Each iteration: error ≈ C * x/T where T is the height of zeros used
# After iteration: we get a new x, recompute, iterate
#
# The question: does this CONVERGE? Or oscillate?

print("\n--- APPROACH 7: Self-consistent iteration with zeta zeros ---")

# First few zeta zeros (imaginary parts)
ZETA_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609809, 65.112544048081607,
    67.079810529494174, 69.546401711173980, 72.067157674481907,
    75.704690699083933, 77.144840068874805,
]

def pi_explicit(x, K=20):
    """Approximate π(x) using explicit formula with K zeta zeros."""
    if x < 2:
        return 0
    result = li(x) - math.log(2)

    # Contribution from zeros
    lnx = math.log(x)
    sqrtx = math.sqrt(x)
    for k in range(min(K, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[k]
        # li(x^ρ) where ρ = 1/2 + iγ
        # Real part of li(x^{1/2+iγ}) ≈ 2*sqrt(x) * cos(γ*ln(x) - angle) / |ρ|*ln(x)
        # More precisely: x^ρ = sqrt(x) * e^{iγ ln x}
        # li(x^ρ) ≈ x^ρ / (ρ * ln x) for large x
        # Real part ≈ sqrtx * cos(gamma*lnx) / (0.5*lnx) ... approximate

        # Better: use the Ei-based formula
        # li(x^ρ) = Ei(ρ * ln x)
        # For complex arg: Ei(a+bi) is complex
        # We need Re[li(x^{1/2+iγ}) + li(x^{1/2-iγ})] = 2 Re[li(x^{1/2+iγ})]

        a = 0.5 * lnx  # real part of ρ*ln(x)
        b = gamma * lnx  # imag part

        # Approximate: li(x^ρ) ≈ x^ρ / (ρ ln x) = sqrt(x) e^{iγ lnx} / ((0.5+iγ) lnx)
        # Real part of pair: 2 Re[sqrt(x) e^{iγ lnx} / ((0.5+iγ) lnx)]
        cos_b = math.cos(b)
        sin_b = math.sin(b)
        denom_real = 0.5 * lnx
        denom_imag = gamma * lnx
        denom_sq = denom_real**2 + denom_imag**2

        re_num = sqrtx * cos_b
        im_num = sqrtx * sin_b

        # (re_num + i*im_num) / (denom_real + i*denom_imag)
        re_result = (re_num * denom_real + im_num * denom_imag) / denom_sq

        result -= 2 * re_result  # factor 2 for conjugate pair

    return result

# Test accuracy of explicit formula with different numbers of zeros
print("\nExplicit formula accuracy (n -> π(p(n)) should equal n):")
test_ns = [100, 500, 1000, 5000, 10000]
for n in test_ns:
    pn = primes[n-1]
    for K in [5, 10, 20]:
        pi_approx = pi_explicit(pn, K)
        error = abs(pi_approx - n)
        print(f"  n={n:6d}, K={K:2d}: π̃(p(n))={pi_approx:.2f}, error={error:.2f}")

# Self-consistent iteration
print("\n--- Self-consistent iteration test ---")
def find_prime_selfconsistent(n, K=20, max_iter=50):
    """Try to find p(n) by iterating with explicit formula."""
    x = R_inverse(n)
    history = [x]

    for it in range(max_iter):
        pi_x = pi_explicit(x, K)
        correction = (n - pi_x) * math.log(x)
        x_new = x + correction
        history.append(x_new)

        if abs(correction) < 0.01:
            break
        x = x_new

    return round(x), history

# Test self-consistent iteration
print("\nSelf-consistent iteration results:")
correct = 0
total = 0
errors = []
for n in range(2, 1001):
    result, history = find_prime_selfconsistent(n, K=20)
    actual = primes[n-1]
    if result == actual:
        correct += 1
    else:
        errors.append((n, result, actual, result - actual))
    total += 1

print(f"Correct: {correct}/{total} ({100*correct/total:.1f}%)")
if errors[:10]:
    print("First 10 errors:")
    for n, got, expected, diff in errors[:10]:
        print(f"  n={n}: got {got}, expected {expected}, diff={diff}")

# APPROACH 8: HIGHER-ORDER EXPLICIT FORMULA
# Instead of π(x), use Π(x) = Σ_{p^k ≤ x} 1/k (Riemann's prime power counting function)
# Π(x) = li(x) - Σ_ρ li(x^ρ) - ln 2 + ∫... (EXACT, no error term if all zeros used)
#
# Then π(x) = Σ_{k=1}^{⌊log₂x⌋} μ(k)/k · Π(x^{1/k})
# = Π(x) - Π(x^{1/2})/2 - Π(x^{1/3})/3 - ...
#
# The key: Π(x) with K zeros has error C*x/(T*ln x)
# But Π(x^{1/2}) with K zeros has error C*x^{1/2}/(T*ln(x^{1/2})) -- MUCH SMALLER!
# So the Möbius inversion REDUCES the error!
#
# Can we exploit this? The total error in π(x) is dominated by the Π(x) term.

print("\n" + "=" * 70)
print("APPROACH 8: Double inversion - explicit formula + Möbius")
print("=" * 70)

def Pi_explicit(x, K=20):
    """Riemann's Π(x) using explicit formula with K zeros."""
    if x < 2:
        return 0
    result = li(x) - math.log(2)

    lnx = math.log(x)
    sqrtx = math.sqrt(x)
    for k in range(min(K, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[k]
        a = 0.5 * lnx
        b = gamma * lnx
        cos_b = math.cos(b)
        sin_b = math.sin(b)
        denom_real = 0.5 * lnx
        denom_imag = gamma * lnx
        denom_sq = denom_real**2 + denom_imag**2
        re_num = sqrtx * cos_b
        im_num = sqrtx * sin_b
        re_result = (re_num * denom_real + im_num * denom_imag) / denom_sq
        result -= 2 * re_result

    return result

def pi_from_Pi(x, K=20):
    """π(x) via Möbius inversion of Π(x)."""
    result = 0
    max_k = int(math.log2(x)) + 1
    for k in range(1, max_k + 1):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        xk = x ** (1.0 / k)
        if xk < 2:
            break
        result += mu_k / k * Pi_explicit(xk, K)
    return result

print("\nπ(x) via Möbius inversion of explicit Π(x):")
for n in [100, 500, 1000, 5000, 10000]:
    pn = primes[n-1]
    for K in [5, 10, 20]:
        pi_val = pi_from_Pi(pn, K)
        error = abs(pi_val - n)
        print(f"  n={n:6d}, K={K:2d}: π(p(n))={pi_val:.2f}, error={error:.2f}")

print("\n--- Summary of exploration ---")
print("This script tests multiple novel approaches.")
print("Key question: can ANY approach avoid the O(sqrt(x)) zeros barrier?")
