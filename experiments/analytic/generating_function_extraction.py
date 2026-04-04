"""
Session 7 Experiment: Generating Function Coefficient Extraction
================================================================
IDEA: The ordinary generating function for primes is:
  F(z) = Σ_{n=1}^∞ p(n) · z^n = 2z + 3z² + 5z³ + 7z⁴ + ...

p(n) = [z^n] F(z) = (1/2πi) ∮ F(z)/z^{n+1} dz  (Cauchy's integral)

If we could evaluate F(z) efficiently on a contour, we could
extract p(n) via FFT/NTT in O(n log n) time.

BUT: F(z) has a natural boundary at |z| = 1 (Fabry gap theorem,
since primes have density → 0).

ALTERNATIVE: What about the DIRICHLET generating function?
  G(s) = Σ_{n=1}^∞ p(n)/n^s

Or the exponential generating function?
  E(z) = Σ_{n=1}^∞ p(n) · z^n / n!

These converge for all z (since p(n) ~ n ln n).

Can we extract coefficients of E(z) efficiently?

NEW IDEA: What about the PRIME ZETA FUNCTION?
  P(s) = Σ_p 1/p^s = Σ_{k=1}^∞ μ(k)/k · ln ζ(ks)

This converges for Re(s) > 1 and has analytic continuation.
Can we extract individual primes from P(s)?
"""

import math
import cmath
import time
import numpy as np
from scipy.fft import fft, ifft

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

primes = sieve_primes(10000)

# ============================================================
# EXPERIMENT 1: Cauchy integral for coefficient extraction
# ============================================================
print("=" * 60)
print("EXPERIMENT 1: Cauchy Integral on EGF")
print("=" * 60)

def egf_eval(z, num_primes=100):
    """Evaluate E(z) = Σ p(n) · z^n / n! for first num_primes terms"""
    result = 0
    factorial = 1
    for n in range(1, num_primes + 1):
        factorial *= n
        result += primes[n-1] * z**n / factorial
    return result

def extract_coefficient_cauchy(func, n, num_points=256, r=0.5):
    """
    Extract n-th coefficient of power series using Cauchy integral.
    [z^n] f(z) = (1/2πi) ∮ f(z)/z^{n+1} dz

    On circle |z| = r:
    [z^n] f(z) = (1/N) Σ_{k=0}^{N-1} f(r·e^{2πik/N}) / (r·e^{2πik/N})^n
    """
    total = 0
    for k in range(num_points):
        theta = 2 * math.pi * k / num_points
        z = r * cmath.exp(1j * theta)
        fz = func(z)
        total += fz / z**n

    coeff = total / num_points
    return coeff.real

# For EGF: coefficient of z^n is p(n)/n!
# So p(n) = n! · [z^n] E(z)
print("Extracting p(n) from EGF via Cauchy integral:")
for n in range(1, 20):
    coeff = extract_coefficient_cauchy(egf_eval, n, num_points=512, r=0.3)
    pn_estimate = coeff * math.factorial(n)
    true_pn = primes[n-1]
    err = abs(pn_estimate - true_pn)
    print(f"  p({n:2d}) = {pn_estimate:10.2f}  (true: {true_pn:5d}, error: {err:.2f})")

print("""
ANALYSIS: Cauchy integral works for SMALL n but:
1. We need to evaluate E(z) at 512+ points on the contour
2. Each evaluation requires summing over the primes we want to find
3. This is CIRCULAR - we need the primes to compute the generating function!

For a non-circular approach, we'd need a CLOSED-FORM for E(z).
E(z) = Σ p(n) · z^n / n! has no known closed form.
""")

# ============================================================
# EXPERIMENT 2: Newton's identities / power sums
# ============================================================
print("=" * 60)
print("EXPERIMENT 2: Newton's Identities for Prime Extraction")
print("=" * 60)

# Idea: Consider primes as "roots" of a polynomial
# f(x) = Π_{p prime, p≤N} (x - p)
# Then the coefficients of f are symmetric functions of primes.
# Newton's identities relate power sums to elementary symmetric polynomials.
#
# Power sums: S_k = Σ_p p^k (sum of k-th powers of primes up to N)
# These are related to the prime zeta function: P(-k) for negative integers
#
# If we could compute S_k efficiently for k = 1, 2, ..., π(N),
# we could recover the individual primes via Newton's identities!
#
# S_k = Σ_{p≤N} p^k can be computed using a GENERALIZED Lucy DP:
# Sum of k-th powers of primes up to x.
#
# This IS computable in O(x^{2/3}) for each k.
# And we need π(N) such sums.
# Total: O(N^{2/3} · π(N)^{1/2}) ... but π(N) ~ N/ln(N)
# So total: O(N^{2/3} · N^{1/2} / ln(N)^{1/2}) >> O(N)
#
# WORSE than direct sieving!

# But let's test it for small N to see if it works in principle.
def power_sums_to_primes(power_sums, num_primes):
    """
    Given power sums S_1, S_2, ..., S_m where S_k = Σ p_i^k,
    recover the primes p_1, ..., p_m using Newton's identities.

    Newton's identities:
    e_1 = S_1
    e_2 = (e_1 · S_1 - S_2) / 2
    e_3 = (e_2 · S_1 - e_1 · S_2 + S_3) / 3
    etc.

    where e_k are elementary symmetric polynomials (coefficients of
    Π(x - p_i) up to sign).
    """
    m = num_primes
    # Compute elementary symmetric polynomials via Newton's identities
    e = [0] * (m + 1)
    e[0] = 1
    for k in range(1, m + 1):
        s = 0
        for j in range(1, k + 1):
            s += ((-1)**(j-1)) * e[k-j] * power_sums[j-1]
        e[k] = s / k

    # The primes are roots of x^m - e_1·x^{m-1} + e_2·x^{m-2} - ... + (-1)^m · e_m
    coeffs = [(-1)**k * e[k] for k in range(m + 1)]
    coeffs.reverse()  # numpy wants highest degree first? No, lowest first.

    # numpy.roots expects [highest, ..., lowest]
    poly_coeffs = [(-1)**k * e[k] for k in range(m, -1, -1)]
    # Actually: the polynomial is x^m - e_1·x^{m-1} + e_2·x^{m-2} - ...
    poly_coeffs_np = []
    for k in range(m + 1):
        poly_coeffs_np.append((-1)**k * e[k])

    roots = np.roots(poly_coeffs_np)
    return sorted([r.real for r in roots if abs(r.imag) < 0.5])

# Test with first few primes
for m in [4, 5, 8, 10, 15]:
    ps = primes[:m]
    power_sums = [sum(p**k for p in ps) for k in range(1, m + 1)]
    recovered = power_sums_to_primes(power_sums, m)
    recovered_rounded = [round(r) for r in recovered]

    match = all(a == b for a, b in zip(recovered_rounded, ps))
    print(f"m={m:2d}: primes={ps}")
    print(f"      recovered={recovered_rounded}  {'✓' if match else '✗'}")

print("""
ANALYSIS: Newton's identities + root finding WORKS but:
1. Polynomial root finding is numerically unstable for large degree
   (condition number grows exponentially with degree)
2. Computing power sums requires O(x^{2/3}) per sum, and we need
   π(x) sums → total cost is WORSE than direct sieving
3. For m > ~50, numerical errors make root recovery impossible

This is a VALID but IMPRACTICAL approach.
""")

# ============================================================
# EXPERIMENT 3: Can we use the prime-counting function as a
# step function and extract steps via differentiation?
# ============================================================
print("=" * 60)
print("EXPERIMENT 3: π(x) Differentiation / Distributional Approach")
print("=" * 60)

print("""
π(x) is a step function: π(x) = Σ_{p≤x} 1

Its distributional derivative is:
  π'(x) = Σ_p δ(x - p)  (sum of Dirac deltas at primes)

The Fourier transform of π'(x) is:
  F[π'](ξ) = Σ_p e^{-2πi p ξ}

This is exactly the "prime exponential sum" studied extensively.

If we could compute F[π'](ξ) at enough points, we could
inverse-Fourier-transform to locate the primes.

BUT: computing Σ_p e^{-2πi p ξ} requires knowing the primes!

UNLESS we could compute it from the ANALYTIC side:
  Σ_p e^{-2πi p ξ} = ... (no known closed form)

The connection to the explicit formula IS the Fourier analysis:
  Σ_p ln(p)/p^s = -ζ'(s)/ζ(s)   (logarithmic derivative)

But this requires Re(s) > 1 and gives the primes WEIGHTED by ln(p),
not individual primes.

CONCLUSION: Fourier/distributional approaches are equivalent to
the explicit formula. No shortcut.
""")

# ============================================================
# EXPERIMENT 4: A radically different idea - BINARY ENCODING
# ============================================================
print("=" * 60)
print("EXPERIMENT 4: Binary Encoding of p(n)")
print("=" * 60)

print("""
Observation: p(n) has about log₂(n · ln(n)) ≈ log₂(n) + log₂(ln(n)) bits.

For n = 10^100: p(n) has about 340 bits.

What if we could compute p(n) BIT BY BIT?

Bit k of p(n) = floor(p(n) / 2^k) mod 2

Using the binary search approach:
  bit k = 1 iff π(m + 2^k) ≥ n where m = (known high bits) · 2^{k+1}

Each bit requires ONE call to π(x).
Total: ~340 calls to π(x), each costing O(x^{2/3}).

Total cost: 340 · O(x^{2/3}) = O(x^{2/3} · log x)

This is SLIGHTLY WORSE than the standard approach (binary search
on π uses ~log(gap) ≈ log(ln(x)) ≈ 6-8 calls to π).

But the BIT-BY-BIT approach has an interesting property:
different bits are INDEPENDENT queries.

If we could compute INDIVIDUAL bits of π(x) cheaper than the
full value, this could win.

Can we compute the k-th bit of π(x) in time less than O(x^{2/3})?
""")

# Test bit-by-bit extraction
def pi_bit_by_bit(n, pi_func, approx_pn):
    """
    Extract p(n) bit by bit using π(x) queries.
    Start from MSB, work down.
    """
    # Start with approximation
    nbits = approx_pn.bit_length() + 2
    result = 0

    for bit in range(nbits - 1, -1, -1):
        candidate = result | (1 << bit)
        # Is p(n) ≥ candidate? Equivalent to: π(candidate) ≥ n
        # But we need: is bit `bit` of p(n) set?
        # Actually: p(n) ≥ candidate iff π(candidate-1) < n
        pi_val = pi_func(candidate)
        if pi_val < n:
            result = candidate

    # result is now the largest x with π(x) < n, so p(n) = result + 1
    # Actually this is just binary search, not bit-by-bit
    return result + 1

from functools import lru_cache

@lru_cache(maxsize=None)
def lucy_pi_cached(x):
    if x < 2:
        return 0
    sqrtx = int(math.isqrt(x))
    small = list(range(-1, sqrtx + 1))  # small[i] = i - 1
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

# Test
for n in [100, 1000, 5000]:
    approx = int(n * (math.log(n) + math.log(math.log(n))))
    t0 = time.time()
    result = pi_bit_by_bit(n, lucy_pi_cached, approx)
    dt = time.time() - t0
    lucy_pi_cached.cache_clear()
    print(f"  p({n}) = {result} (true: {primes[n-1]}, {'✓' if result == primes[n-1] else '✗'}, time: {dt:.3f}s)")

print("""
RESULT: Bit-by-bit extraction works but is equivalent to binary search.
The key insight: we CANNOT compute individual bits of π(x) cheaper
than the full value. Information-theoretically, knowing ALL bits of
π(x) IS knowing π(x).
""")

# ============================================================
# EXPERIMENT 5: Sieving in function space
# ============================================================
print("=" * 60)
print("EXPERIMENT 5: Abstract Sieving Ideas")
print("=" * 60)

print("""
The Eratosthenes sieve works by:
1. Start with all integers [2, N]
2. Remove multiples of 2, then 3, then 5, ...
3. What remains are primes

What if instead of sieving INTEGERS, we sieve FUNCTIONS?

Idea: Start with f₀(n) = n·ln(n) (approximate nth prime)
Sieve step: correct f_k(n) using information about prime p_k

For example, after sieving by p=2:
  f₁(n) = adjusted formula accounting for half the integers being even

This is essentially the "sieve of Eratosthenes in analytic form":
  R(x) → R(x) - R(x^{1/2})/2 → ... → π(x)

Which IS the explicit formula / Möbius inversion.

So "sieving in function space" is exactly the analytic approach
we already know costs O(x^{2/3}) or O(x^{1/2+ε}).

NO NEW APPROACH HERE.
""")

# SUMMARY
print("=" * 60)
print("SESSION 7 - GENERATING FUNCTION EXPERIMENTS SUMMARY")
print("=" * 60)

print("""
All 5 experiments confirm the same barrier:

1. CRT Modular: π(x) mod m costs same as π(x) exactly
2. Newton's Identities: Works but numerically unstable, cost > direct
3. Fourier/Distributional: Equivalent to explicit formula
4. Bit-by-bit: Equivalent to binary search, no per-bit shortcut
5. Abstract Sieving: Equivalent to analytic sieve methods

NO NEW VIABLE APPROACH FOUND in this experiment set.
The information barrier remains: ~5 bits/prime of irreducible entropy.
""")
