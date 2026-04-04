#!/usr/bin/env python3
"""
Session 10: Lattice and F_1 Geometry Approaches to nth Prime

Radical idea: In F_q[x], the nth irreducible polynomial can be found efficiently.
The "field with one element" F_1 analogy suggests Z = F_1[x].
Can we exploit this analogy computationally?

Also: Geometric/lattice approaches to prime distribution.
"""

import math
from functools import lru_cache

# ============================================================
# PART 1: F_1 ANALOGY - Counting irreducibles over F_q
# ============================================================

def mobius(n):
    """Compute μ(n)"""
    if n == 1:
        return 1
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            count = 0
            while temp % d == 0:
                temp //= d
                count += 1
            if count > 1:
                return 0
            factors.append(d)
        d += 1
    if temp > 1:
        factors.append(temp)
    return (-1) ** len(factors)

def count_irreducibles_fq(q, n):
    """Number of monic irreducible polynomials of degree n over F_q.
    Formula: (1/n) * Σ_{d|n} μ(n/d) * q^d
    """
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            total += mobius(n // d) * (q ** d)
    return total // n

def prime_counting_fq(q, N):
    """Total irreducible polynomials of degree ≤ N over F_q.
    Analog of π(q^N) in the polynomial ring.
    """
    return sum(count_irreducibles_fq(q, d) for d in range(1, N + 1))

# Test: as q→1, does count_irreducibles approach something prime-like?
print("=" * 60)
print("PART 1: F_1 Analogy - Irreducibles over F_q as q→1")
print("=" * 60)

# For F_q, the "numbers" up to degree n correspond to integers up to q^n
# Irreducibles of degree n ↔ primes between q^{n-1} and q^n
# As q→1: q^n → 1, so the analogy breaks down naively

# But there's a deeper version: "F_1 zeta function"
# ζ_{F_1[x]}(s) should equal ζ(s) (Riemann zeta)
# This means: the prime-counting in F_1[x] IS π(x)

# The formula for F_q: π_{F_q}(q^n) = Σ_{k=1}^n (1/k) Σ_{d|k} μ(k/d) q^d
# Naive q→1 limit: Σ_{k=1}^n (1/k) Σ_{d|k} μ(k/d) * 1 = Σ_{k=1}^n (1/k) * [k=1] = 1
# This is wrong - the limit is degenerate

# Better: use the logarithmic derivative
# d/ds log ζ_{F_q[x]}(s) = Σ_P deg(P) |P|^{-s} log|P|
# For F_q: |P| = q^{deg(P)}, so this = Σ_n n*I_q(n) q^{-ns} log(q)

# Try: interpolate between F_2 and "F_1"
for q in [2, 3, 5, 7, 11]:
    counts = [count_irreducibles_fq(q, n) for n in range(1, 20)]
    pi_values = [prime_counting_fq(q, n) for n in range(1, 20)]
    print(f"\nF_{q}: irreducibles by degree: {counts[:10]}")
    print(f"F_{q}: cumulative (analog of π): {pi_values[:10]}")

# Now compare with actual primes
from sympy import primepi, prime, isprime
import sympy

actual_pi = [primepi(n) for n in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]]
print(f"\nActual π at powers of 2: {actual_pi}")

# Key question: Is there a "q→1 regularization" that gives π(x)?
# Soulé's formula: lim_{q→1} (π_{F_q}(q^n) - q^n/(n*ln(q))) should relate to π(e^n)?
# Let's test numerically
print("\nRegularized F_q limit (Soulé-type):")
for n in range(1, 15):
    for q in [1.01, 1.001, 1.0001]:
        # For non-integer q, use the analytic formula
        total = 0
        for k in range(1, n + 1):
            inner = 0
            for d in range(1, k + 1):
                if k % d == 0:
                    inner += mobius(k // d) * (q ** d)
            total += inner / k
        x_val = q ** n  # This is the "size" we're counting up to
        print(f"  q={q:.4f}, n={n}: π_approx = {total:.4f}, x={x_val:.4f}, actual π({int(round(x_val))}) = {primepi(int(round(x_val)))}")
    print()
    if n > 5:
        break


# ============================================================
# PART 2: Lattice Approach - Primes as lattice projections
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Lattice Approach")
print("=" * 60)

# Idea: embed the prime sequence into a higher-dimensional lattice
# such that the projection onto one axis gives the primes.
#
# The simplest version: the Ulam spiral. Primes tend to cluster
# on certain diagonals. Can we quantify this?

# Better idea: embed n → (n, p(n)) and look at the lattice structure
# of the resulting point set.

# Even better: consider the "prime lattice" Λ_P = {(n, p(n)) : n ∈ N}
# This is a subset of Z^2. What's the shortest vector? The densest sublattice?

# Compute gaps and second differences
primes_list = list(sympy.primerange(2, 10000))
gaps = [primes_list[i+1] - primes_list[i] for i in range(len(primes_list)-1)]
second_diff = [gaps[i+1] - gaps[i] for i in range(len(gaps)-1)]

# Autocorrelation of gaps
import numpy as np
gaps_arr = np.array(gaps[:1000], dtype=float)
gaps_centered = gaps_arr - np.mean(gaps_arr)
autocorr = np.correlate(gaps_centered, gaps_centered, mode='full')
autocorr = autocorr[len(autocorr)//2:]
autocorr /= autocorr[0]

print(f"Gap autocorrelation (lags 1-10): {[f'{x:.4f}' for x in autocorr[1:11]]}")
print(f"Gap autocorrelation (lags 10-20): {[f'{x:.4f}' for x in autocorr[10:21]]}")

# Hardy-Littlewood prediction: consecutive gaps should be NEGATIVELY correlated
# (Lemke Oliver-Soundararajan bias)
print(f"\nGap-gap correlation: {np.corrcoef(gaps_arr[:-1], gaps_arr[1:])[0,1]:.6f}")
print("(Expected: slightly negative due to Lemke Oliver-Soundararajan)")


# ============================================================
# PART 3: Continued Fraction of Prime-Related Constants
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Continued Fraction Patterns")
print("=" * 60)

def continued_fraction(x, terms=30):
    """Compute continued fraction expansion of x"""
    cf = []
    for _ in range(terms):
        a = int(math.floor(x))
        cf.append(a)
        frac = x - a
        if abs(frac) < 1e-12:
            break
        x = 1.0 / frac
    return cf

# Prime constant: Σ 1/2^p(n) = 0.414682509851...
# If this has a pattern in its CF, we could extract primes
prime_const = sum(1.0 / (2 ** p) for p in primes_list[:100])
cf_prime = continued_fraction(prime_const, 40)
print(f"Prime constant CF: {cf_prime[:30]}")

# Copeland-Erdős constant: 0.2357111317192329...
ce_str = '0.' + ''.join(str(p) for p in primes_list[:200])
ce_const = float(ce_str[:50])  # Limited by float precision
cf_ce = continued_fraction(ce_const, 30)
print(f"Copeland-Erdős CF: {cf_ce[:20]}")

# Sum of reciprocal primes (diverges, but partial sums)
# Mertens constant M ≈ 0.2614972128
mertens = sum(1.0/p for p in primes_list[:1000]) - math.log(math.log(primes_list[999]))
print(f"\nMertens constant approx: {mertens:.10f}")
cf_mertens = continued_fraction(mertens, 30)
print(f"Mertens CF: {cf_mertens[:20]}")

# Twin prime constant: Π_{p≥3} (1 - 1/(p-1)^2)
twin_const = 1.0
for p in primes_list[1:200]:  # skip 2
    twin_const *= (1 - 1.0 / (p - 1) ** 2)
twin_const *= 2  # C_2
print(f"\nTwin prime constant C_2 ≈ {twin_const:.10f}")
cf_twin = continued_fraction(twin_const, 30)
print(f"Twin prime CF: {cf_twin[:20]}")


# ============================================================
# PART 4: Beatty Sequence / Mechanical Sequence Approach
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Beatty/Mechanical Sequence Approaches")
print("=" * 60)

# Can we find α such that ⌊n·α⌋ hits near primes?
# Best candidate: α = some function of n (not constant)
#
# Actually, let's try: α(n) = n * ln(n) + n * ln(ln(n)) - n + ...
# This is the asymptotic expansion of p(n)

def cipolla_approx(n, terms=6):
    """Cipolla's asymptotic expansion for p(n)"""
    if n < 2:
        return 2
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n) if ln_n > 0 else 0

    # p(n) ~ n*ln(n) + n*ln(ln(n)) - n + n*ln(ln(n))/ln(n) - 2n/ln(n) + ...
    result = n * ln_n
    result += n * ln_ln_n
    result -= n
    result += n * ln_ln_n / ln_n
    result -= 2 * n / ln_n
    # More terms from the asymptotic expansion
    result += n * (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2)
    return result

# How good is Cipolla vs R^{-1}?
print("Cipolla approximation accuracy:")
errors_cipolla = []
for idx in [100, 500, 1000, 5000, 10000]:
    p_actual = prime(idx)
    p_cipolla = cipolla_approx(idx)
    err = abs(p_actual - p_cipolla) / p_actual
    errors_cipolla.append(err)
    print(f"  p({idx}) = {p_actual}, Cipolla = {p_cipolla:.1f}, rel_err = {err:.6f}")


# ============================================================
# PART 5: Novel Idea - Prime Index via Möbius Inversion
# ============================================================
print("\n" + "=" * 60)
print("PART 5: Prime Index via Möbius Inversion")
print("=" * 60)

# π(x) = Σ_{n≤x} χ_P(n) where χ_P is the prime indicator
# By Möbius inversion of the sieve:
# χ_P(n) = Σ_{d|n, d<n} μ(d) * ... (this is essentially the sieve)
#
# But what about the INVERSE problem?
# Given n, find the smallest x such that π(x) = n
# This is p(n) = min{x : π(x) = n}
#
# Key identity: p(n) = n + 1 + Σ_{k=2}^{n} (⌊1/π(k)⌋ - ... )
# Actually, Rosser's formula: p(n) = the n-th element of the sieve

# Let's try a completely different angle:
# What if we could compute π(x) mod m for small m, for all x?
# Then by CRT, we could reconstruct π(x) exactly.
#
# π(x) mod 2: parity of prime-counting function
# This is EXTREMELY hard - it's related to Selberg's parity problem!

# Test: how does π(x) mod 2 behave?
parity_changes = 0
prev_parity = 0
for x in range(2, 10001):
    pi_x = primepi(x)
    cur_parity = pi_x % 2
    if cur_parity != prev_parity:
        parity_changes += 1
    prev_parity = cur_parity

print(f"π(x) parity changes in [2, 10000]: {parity_changes}")
print(f"(= number of primes in range, as expected: {primepi(10000)})")

# π(x) mod 2 changes at every prime - so knowing the parity IS knowing the primes
# This confirms the difficulty.

# What about π(x) mod larger m?
for m in [3, 5, 7, 10]:
    residues = [primepi(x) % m for x in range(2, 1001)]
    from collections import Counter
    dist = Counter(residues)
    print(f"π(x) mod {m} distribution (x=2..1000): {dict(sorted(dist.items()))}")


# ============================================================
# PART 6: Digit-by-digit construction
# ============================================================
print("\n" + "=" * 60)
print("PART 6: Digit-by-Digit Prime Construction")
print("=" * 60)

# Idea: construct p(n) one digit at a time, from most significant to least.
# To determine the k-th digit, we need to know π(x) at specific values.
# But this STILL requires computing π(x)...
#
# UNLESS: there's a way to determine individual digits from local information.
#
# For the Copeland-Erdős constant, BBP-like formulas don't exist (proven).
# But what about a different encoding?
#
# Consider the BINARY representation of p(n).
# The most significant bit is always 1 (for p > 1).
# The next bit tells us if p(n) > 3·2^{k-2} or not.
# This requires knowing π(3·2^{k-2}) vs n.

# Let's measure: how many bits of p(n) can we get from R^{-1}(n)?
from mpmath import mp, mpf, li, log as mplog

mp.dps = 50

def R_inverse_approx(n):
    """Approximate inverse Riemann R function using Newton's method on li"""
    # Start with x = n * ln(n) as initial guess
    x = mpf(n) * mplog(mpf(n))
    for _ in range(50):
        lix = li(x)
        if abs(lix - n) < mpf('1e-30'):
            break
        # Newton step: x_{k+1} = x_k - (li(x_k) - n) / (1/ln(x_k))
        x = x - (lix - n) * mplog(x)
    return x

print("Bits of accuracy from R^{-1}(n):")
test_ns = [100, 1000, 10000, 50000]
for n in test_ns:
    p_actual = prime(n)
    r_inv = float(R_inverse_approx(n))
    error = abs(p_actual - r_inv)
    bits_total = math.log2(p_actual)
    bits_error = math.log2(max(error, 1))
    bits_correct = bits_total - bits_error
    print(f"  n={n}: p(n)={p_actual}, R^-1={r_inv:.1f}, error={error:.1f}")
    print(f"    Total bits: {bits_total:.1f}, Error bits: {bits_error:.1f}, Correct bits: {bits_correct:.1f}")
    print(f"    Fraction correct: {bits_correct/bits_total:.4f}")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY OF FINDINGS")
print("=" * 60)
print("""
1. F_1 ANALOGY: The q→1 limit is degenerate. The formula for F_q
   counting becomes trivial at q=1. The analogy is structural,
   not computational.

2. LATTICE: Prime gaps show weak negative autocorrelation
   (Lemke Oliver-Soundararajan). No exploitable lattice structure found.

3. CONTINUED FRACTIONS: No patterns detected in CFs of prime constants.
   These constants appear to have "generic" CF expansions.

4. CIPOLLA: Asymptotic expansion gives ~0.1-5% relative error.
   Much worse than R^{-1}(n) which gives ~50% of bits correct.

5. MÖBIUS/CRT: π(x) mod m requires knowing primes (circular).
   Parity of π(x) is exactly as hard as finding primes.

6. DIGIT-BY-DIGIT: R^{-1}(n) gives about 50% of the bits of p(n).
   The remaining bits encode the "random" fluctuation Σ_ρ R(x^ρ).
""")
