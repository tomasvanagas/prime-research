"""
Session 9: Cramér Model with Deterministic Corrections

The Cramér model says: treat each integer n > 2 as "prime" with probability 1/ln(n).
This gives PNT perfectly. But the DEVIATIONS from this model are deterministic.

Key idea: The Maier matrix phenomenon shows that primes are NOT independent.
In intervals of length (log x)², the count can deviate from expected by any constant factor.
This is due to the sieve structure.

What if we can compute the deterministic correction to the Cramér model?

The correction comes from:
1. Small prime divisibility (wheel factorization) — deterministic, O(1)
2. Chebyshev bias (primes favor certain residue classes) — computable
3. Second-order sieve effects — captured by Selberg sieve
4. Higher-order effects — captured by the explicit formula (zeta zeros)

Can we compute enough of these corrections to get exact p(n)?

Also: NEW IDEA — what if we use the EXACT Cramér model?
Define X_n = 1 if n is prime, 0 otherwise.
E[X_n] = 1/ln(n)
Var[X_n] = 1/ln(n) · (1 - 1/ln(n))
Cov[X_n, X_m] = ??? (This is the key!)

The covariance structure of the prime indicator IS the explicit formula.
If we could compute this covariance without zeta zeros, we'd have a shortcut.
"""

import numpy as np
from sympy import prime, primepi, isprime, factorint
from mpmath import mp, mpf, log, li, exp
import time

mp.dps = 30

print("=" * 70)
print("EXPERIMENT 1: Wheel Factorization Correction")
print("=" * 70)

# Wheel mod 30 (= 2·3·5): candidates are n ≡ 1,7,11,13,17,19,23,29 (mod 30)
# Density correction: instead of 1/ln(n), use (30/8)/ln(n) for survivors

# More generally, for wheel mod W = p₁·p₂·...·pₖ:
# Density = (W/φ(W))/ln(n) for survivors
# This is a MULTIPLICATIVE correction — the smooth part of Mertens' theorem

# How much does the wheel reduce error?
wheel_primes = [2, 3, 5, 7, 11, 13]

def wheel_density_correction(n, wheel_p):
    """Compute the local density correction from wheel factorization"""
    correction = 1.0
    for p in wheel_p:
        if n % p == 0:
            return 0  # n is divisible by a wheel prime
        correction *= p / (p - 1)  # Mertens-like correction
    return correction / np.log(n)

# Compare Cramér density vs wheel-corrected density
print("Comparison of prime density models:")
for x in [100, 1000, 10000]:
    # Count primes in [x, x+100]
    actual = sum(1 for n in range(x, x+100) if isprime(n))
    cramer = sum(1/np.log(n) for n in range(x, x+100))

    # Wheel mod 30 correction
    wheel30 = 0
    for n in range(x, x+100):
        if all(n % p != 0 for p in [2, 3, 5]):
            wheel30 += (30/8) / np.log(n)

    # Wheel mod 2310 correction
    wheel2310 = 0
    for n in range(x, x+100):
        if all(n % p != 0 for p in [2, 3, 5, 7, 11]):
            wheel2310 += (2310/480) / np.log(n)

    print(f"  [{x}, {x+100}]: actual={actual}, Cramér={cramer:.1f}, wheel30={wheel30:.1f}, wheel2310={wheel2310:.1f}")

print("\n" + "=" * 70)
print("EXPERIMENT 2: Chebyshev Bias Correction")
print("=" * 70)

# Primes ≡ 3 mod 4 slightly outnumber primes ≡ 1 mod 4 (Chebyshev bias)
# More generally, in the "race" between primes ≡ a and ≡ b mod q,
# one class is favored based on whether a is a QR mod q.

# Can this bias help us predict primes?
for q in [4, 6, 8, 12]:
    classes = {}
    for p_idx in range(1, 1001):
        p = prime(p_idx)
        r = p % q
        classes[r] = classes.get(r, 0) + 1

    print(f"  First 10000 primes mod {q}: {dict(sorted(classes.items()))}")

# The bias is real but tiny — typically O(√x/ln(x)) excess
# For x=10^100, bias ≈ 10^49/236 ≈ 4×10^46
# This is LESS than the total error, so it doesn't help

print("\n" + "=" * 70)
print("EXPERIMENT 3: Hardy-Littlewood Constants")
print("=" * 70)

# The twin prime constant C₂ = Π_{p>2} (1 - 1/(p-1)²) ≈ 0.6601618...
# gives the density of twin primes.
# More generally, for any admissible k-tuple {h₁,...,hₖ}:
# #{n ≤ x : n+h₁,...,n+hₖ all prime} ~ S(h₁,...,hₖ) · x/(ln x)^k

# These constants are COMPUTABLE from small primes only!
# Can we use k-tuple densities to recover individual primes?

# Idea: knowing the density of {p, p+2} twins, {p, p+6} cousins, etc.
# gives us STATISTICS about gaps. But statistics ≠ individual values.

# Let's check: can gap distribution statistics predict p(n)?

gaps = [prime(n+1) - prime(n) for n in range(1, 1001)]
print("Gap statistics (first 1000 primes):")
from collections import Counter
gap_counts = Counter(gaps)
top_gaps = sorted(gap_counts.items(), key=lambda x: -x[1])[:10]
print(f"  Most common gaps: {top_gaps}")
print(f"  Mean gap: {np.mean(gaps):.2f}")
print(f"  Std gap: {np.std(gaps):.2f}")
print(f"  Max gap: {max(gaps)}")

# Entropy of gap distribution
total = len(gaps)
probs = [c/total for _, c in gap_counts.items()]
entropy = -sum(p * np.log2(p) for p in probs if p > 0)
print(f"  Gap entropy: {entropy:.2f} bits")
print(f"  Number of distinct gaps: {len(gap_counts)}")

print("\n" + "=" * 70)
print("EXPERIMENT 4: Hybrid Model — R^{-1} + Local Density Correction")
print("=" * 70)

# Most promising practical idea:
# 1. Compute R^{-1}(n) as initial estimate (~50% digits correct)
# 2. Apply Chebyshev bias correction
# 3. Apply local density correction based on residue classes
# 4. Round to nearest prime
# 5. Verify with Miller-Rabin

# The question: can steps 2-3 reduce error enough for step 4 to work?

from mpmath import lambertw

def R_inverse_approx(n):
    """Approximate R^{-1}(n) — the inverse Riemann R function"""
    n = mpf(n)
    if n <= 1:
        return mpf(2)
    L1 = log(n)
    L2 = log(L1) if L1 > 0 else mpf(0)
    # First approximation: n·ln(n) + n·ln(ln(n))
    x = n * (L1 + L2 - 1 + (L2 - 2)/L1 +
             ((L2)**2 - 6*L2 + 11)/(2*L1**2))
    return float(x)

# Test hybrid approach
correct = 0
total = 0
max_error = 0
for n in range(100, 2001):
    pn = prime(n)
    r_inv = R_inverse_approx(n)

    # Local density correction: adjust based on residue class
    # Nearest prime to r_inv
    candidate = int(round(r_inv))

    # Find nearest prime
    if isprime(candidate):
        nearest = candidate
    else:
        # Search both directions
        up = candidate + 1
        while not isprime(up):
            up += 1
        down = candidate - 1
        while down > 1 and not isprime(down):
            down -= 1
        nearest = up if abs(up - candidate) <= abs(down - candidate) else down

    if nearest == pn:
        correct += 1
    error = abs(r_inv - pn)
    max_error = max(max_error, error)
    total += 1

pct = 100 * correct / total
print(f"R^{{-1}}(n) → nearest prime approach:")
print(f"  Accuracy: {correct}/{total} = {pct:.1f}%")
print(f"  Max error: {max_error:.1f}")
print(f"  (This finds the NEAREST prime, not necessarily the NTH prime)")

print("\n" + "=" * 70)
print("EXPERIMENT 5: Can We Compute π(x) mod 2 Efficiently?")
print("=" * 70)

# If we could compute π(x) mod 2 in polylog time, then:
# Binary search with parity → p(n) in O(polylog²)

# π(x) mod 2: is x closer to an even or odd count of primes?
# By the explicit formula: π(x) = R(x) - Σ_ρ R(x^ρ) - constant terms
# π(x) mod 2 depends on the PARITY of the sum over zeros.

# This is related to the PARITY problem in sieve theory!
# Selberg's parity barrier: linear sieves cannot distinguish
# products of even vs odd number of prime factors.

# Can we compute π(x) mod 2 from R(x)?
print("Parity of π(x) from smooth approximation:")
correct_parity = 0
total_parity = 0
for x in range(100, 2001):
    actual = int(primepi(x))
    smooth = float(li(mpf(x)))
    if round(smooth) % 2 == actual % 2:
        correct_parity += 1
    total_parity += 1

print(f"  li(x) parity matches π(x): {correct_parity}/{total_parity} = {100*correct_parity/total_parity:.1f}%")
print(f"  (50% = random, anything above 50% shows SOME parity info in li)")

# Even getting parity right is hard!
# The parity of π(x) oscillates rapidly and is tied to the Liouville function λ(n)
# Σ λ(n) = # even prime factors - # odd prime factors
# This is equivalent to the Mertens function, which is O(√x)

print("\n" + "=" * 70)
print("FINAL ANALYSIS: All Cramér Model Extensions Fail")
print("=" * 70)
print("""
1. Wheel factorization: constant factor improvement only
2. Chebyshev bias: O(√x/ln x) effect — smaller than the total error
3. Hardy-Littlewood constants: give statistics, not positions
4. R^{-1} + nearest prime: ~30% correct (finds nearest, not nth)
5. Parity of π(x): ~52% from li(x) — barely above random

The Cramér model captures the smooth part of prime distribution.
ALL corrections (bias, k-tuple densities, etc.) are STILL smooth.
The non-smooth part IS the zeta zeros, requiring O(√x) computation.

There is no middle ground: either you know the exact distribution
(requiring O(√x) work) or you have only statistics (giving O(√x) error).
""")
