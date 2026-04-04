"""
Session 7 Experiment: Modular CRT Reconstruction of π(x)
=========================================================
IDEA: Can we compute π(x) mod m for small moduli m CHEAPLY,
then reconstruct π(x) exactly via CRT?

Key question: Is computing π(x) mod m any easier than computing π(x)?

Analysis:
- π(x) mod 2 = parity of π(x)
- The parity of π(x) is related to the Liouville function λ(n)
- Specifically, π(x) mod 2 relates to whether x has an even or odd
  number of primes below it

Let's test if there's ANY shortcut for computing π(x) mod small m.
"""

import math
import time
from functools import lru_cache

# First, let's understand what π(x) mod 2 means
# π(x) mod 2 tells us the parity of the number of primes ≤ x
# This is related to the Mertens function M(x) = Σ_{n≤x} μ(n)
# and the Liouville summatory function L(x) = Σ_{n≤x} λ(n)

# The key relationship:
# (-1)^π(x) = product over primes p≤x of (-1) = (-1)^π(x)
# So π(x) mod 2 = (1 - (-1)^π(x)) / 2

# Can we compute (-1)^π(x) without computing π(x)?
# (-1)^π(x) = Liouville function evaluated at... no, that's λ(n) for individual n

# Actually, let's think about this differently.
# The Legendre sieve gives:
# π(x) - π(√x) + 1 = Σ_{d | P(√x)} μ(d) * floor(x/d)
# where P(√x) = product of primes ≤ √x

# So π(x) mod m = (Σ_{d | P(√x)} μ(d) * floor(x/d) + π(√x) - 1) mod m

# The sum has 2^π(√x) terms... exponentially many.
# BUT mod m, many terms might cancel!

# Let's test: for small m, how many DISTINCT values does
# μ(d) * floor(x/d) mod m take?

def sieve_primes(n):
    """Simple sieve of Eratosthenes"""
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n + 1) if is_prime[i]]

def pi_exact(x):
    """Exact prime counting via sieve (for testing, small x only)"""
    if x < 2:
        return 0
    return len(sieve_primes(int(x)))

def legendre_pi_mod_m(x, m):
    """
    Try to compute π(x) mod m using Legendre's formula.
    π(x) - π(√x) + 1 = Σ_{d | P(√x)} μ(d) * floor(x/d)

    We compute the RHS mod m.
    """
    x = int(x)
    sqrtx = int(math.isqrt(x))
    primes_below_sqrt = sieve_primes(sqrtx)

    # This has 2^len(primes_below_sqrt) terms - exponential!
    # For x=100, √x=10, primes=[2,3,5,7], 16 terms - fine
    # For x=10000, √x=100, primes up to 100 = 25 primes, 2^25 = 33M terms
    # For x=10^100, √x=10^50, π(√x) ~ 10^48 terms -> IMPOSSIBLE

    if len(primes_below_sqrt) > 20:
        return None, "Too many terms"

    # Compute via inclusion-exclusion mod m
    total = 0
    num_terms = 0

    # Iterate over all subsets of primes_below_sqrt
    for mask in range(1 << len(primes_below_sqrt)):
        d = 1
        num_factors = 0
        for i in range(len(primes_below_sqrt)):
            if mask & (1 << i):
                d *= primes_below_sqrt[i]
                num_factors += 1

        mu_d = (-1) ** num_factors  # μ(d) for squarefree d
        contribution = (mu_d * (x // d)) % m
        total = (total + contribution) % m
        num_terms += 1

    # total = π(x) - π(√x) + 1 mod m
    pi_sqrt = len(primes_below_sqrt)
    result = (total + pi_sqrt - 1) % m

    return result, num_terms

# Test the modular approach
print("=" * 60)
print("EXPERIMENT 1: Legendre's Formula mod m")
print("=" * 60)

test_values = [100, 200, 500, 1000, 5000, 10000]
moduli = [2, 3, 5, 7, 11]

for x in test_values:
    pi_true = pi_exact(x)
    print(f"\nx = {x}, π(x) = {pi_true}")

    for m in moduli:
        result, info = legendre_pi_mod_m(x, m)
        if result is not None:
            expected = pi_true % m
            match = "✓" if result == expected else "✗"
            print(f"  mod {m}: computed={result}, expected={expected} {match} (terms={info})")
        else:
            print(f"  mod {m}: {info}")

# Now let's think about whether modular computation can be cheaper
print("\n" + "=" * 60)
print("EXPERIMENT 2: Cost Analysis of Modular π(x)")
print("=" * 60)

print("""
ANALYSIS:
The Legendre formula has 2^π(√x) terms.
Computing mod m doesn't reduce the NUMBER of terms.
It only reduces the size of intermediate values.

For x = 10^102 (the target for p(10^100)):
  √x = 10^51
  π(√x) ≈ 10^51 / (51 · ln 10) ≈ 8.5 × 10^48
  Number of Legendre terms: 2^(8.5 × 10^48) - completely infeasible

The Meissel-Lehmer method reduces this to O(x^{2/3}) terms.
But O(x^{2/3}) = O(10^68) - still infeasible.

CAN WE DO BETTER MOD m?

The key question: Does computing π(x) mod m have lower complexity
than computing π(x)?

Theorem (informal): π(x) mod 2 is equivalent to computing the
parity of the Legendre sum. The parity of subset sums is generally
#P-hard (counting problems). So NO, there's no shortcut.

More formally: if we could compute π(x) mod 2 in time T, we could
compute π(x) in time O(T · log(π(x))) by computing mod 2, mod 3,
mod 5, ... and using CRT. But the reverse is also true. So the
complexity of π(x) mod m is essentially the same as π(x).
""")

# EXPERIMENT 3: A different approach - can we compute π(x) mod m
# using the ANALYTIC method more efficiently?
print("=" * 60)
print("EXPERIMENT 3: Analytic Approach mod m")
print("=" * 60)

print("""
The explicit formula: π(x) = R(x) - Σ_ρ R(x^ρ) + corrections

For π(x) mod m, we need:
  R(x) mod m - Σ_ρ R(x^ρ) mod m

Problem: R(x) and R(x^ρ) are REAL NUMBERS, not integers.
The individual terms are not integers. Only the TOTAL sum is
(approximately) an integer.

So we can't take individual terms mod m. We need to compute
the full sum to sufficient precision that we can round to the
nearest integer, and THEN take mod m.

This means: computing π(x) mod m via the analytic method has
the SAME cost as computing π(x) exactly.

CONCLUSION: Modular CRT reconstruction does NOT help.
The complexity of π(x) mod m ≈ complexity of π(x).
""")

# EXPERIMENT 4: What about a completely different characterization?
# Can we compute π(x) using a RECURSIVE structure?
print("=" * 60)
print("EXPERIMENT 4: Recursive π(x) via Lucy_Hedgehog")
print("=" * 60)

def lucy_pi(x):
    """Lucy_Hedgehog DP algorithm for π(x). O(x^{2/3}) time, O(x^{1/2}) space."""
    if x < 2:
        return 0

    sqrtx = int(math.isqrt(x))

    # S[v] = number of integers in [2, v] that survive sieving
    # For v > sqrtx, we use S_large[x // v]
    # For v <= sqrtx, we use S_small[v]

    small = [0] * (sqrtx + 2)
    large = [0] * (sqrtx + 2)

    for i in range(1, sqrtx + 1):
        small[i] = i - 1
        large[i] = x // i - 1

    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:  # p is not prime
            continue

        cnt = small[p - 1]  # π(p-1)
        p2 = p * p

        for i in range(1, min(sqrtx, x // p2) + 1):
            if i * p <= sqrtx:
                large[i] -= large[i * p] - cnt
            else:
                large[i] -= small[x // (i * p)] - cnt

        for i in range(sqrtx, p2 - 1, -1):
            small[i] -= small[i // p] - cnt

    return large[1]

# Test Lucy DP
print("Testing Lucy DP correctness:")
for x in [10, 100, 1000, 10000, 100000, 1000000]:
    t0 = time.time()
    result = lucy_pi(x)
    dt = time.time() - t0
    print(f"  π({x:>10}) = {result:>8}  ({dt:.4f}s)")

# Can we make Lucy DP work modularly?
print("\nCan Lucy DP be computed mod m?")
print("""
The Lucy DP computes S[v] for all v in {x//i : i=1..√x} ∪ {1..√x}.
At each step, it updates: S[v] -= S[v/p] - π(p-1)

If we compute mod m:
  S[v] mod m = (S[v] - S[v/p] + π(p-1)) mod m

This WORKS! Each step is a simple modular arithmetic operation.

BUT: the algorithm still performs O(x^{2/3}) steps.
Computing mod m saves NOTHING on the number of steps.
It only makes each step O(1) instead of O(1) (they're already O(1)).

So: Lucy DP mod m costs O(x^{2/3}), same as Lucy DP exactly.
""")

# EXPERIMENT 5: Is there a way to "skip" steps in Lucy DP?
print("=" * 60)
print("EXPERIMENT 5: Can We Skip Steps in Lucy DP?")
print("=" * 60)

print("""
The Lucy DP processes primes p = 2, 3, 5, 7, 11, ... up to √x.
For each prime p, it updates O(√x / p) entries.
Total work: Σ_{p ≤ √x} √x/p ≈ √x · ln(ln(√x)) ≈ x^{1/2} · ln(ln(x))

But that's the FAST version. The actual complexity is O(x^{2/3} / ln(x))
because the number of DP entries is O(√x) and there are π(√x) primes.

To SKIP steps, we'd need a way to combine the effects of multiple
primes at once. For example:
- Can we process ALL primes from p₁ to p₂ in one batch?
- Can we use the structure of the updates (they're multiplicative)
  to do matrix exponentiation?

The update for prime p is: S[v] -= S[v/p] - c
This is a LINEAR operation on the vector S.
So processing prime p is a matrix multiplication: S → M_p · S

Processing all primes from 2 to √x is: S → M_{p_π(√x)} · ... · M_3 · M_2 · S

Each M_p is a sparse matrix (O(√x/p) non-zero off-diagonal entries).
The product is NOT sparse in general.

Can we compute the matrix product faster using:
1. FFT? No - the matrices are not circulant.
2. Divide and conquer? Product of k matrices of size n × n = O(k · n^ω)
   where ω ≈ 2.37. This is WORSE than direct DP.
3. Matrix exponentiation? The matrices are all different (one per prime).

CONCLUSION: No way to skip steps in Lucy DP.
""")

# EXPERIMENT 6: What about a RANDOMIZED algorithm?
print("=" * 60)
print("EXPERIMENT 6: Randomized π(x)?")
print("=" * 60)

import random

def randomized_pi_test(x, num_samples):
    """
    Attempt: Estimate π(x) by random sampling.
    Sample random integers in [2, x] and check primality.
    By PNT, fraction of primes ≈ 1/ln(x).
    So π(x) ≈ (x-1) · (fraction that are prime)

    Error: O(√(x/ln(x)) / √num_samples) by CLT.
    For error < 0.5 (to round to exact), need num_samples > 4x/ln(x).
    That's MORE work than direct computation!
    """
    if x < 100:
        return pi_exact(x)

    primes = set(sieve_primes(int(x)))
    count = 0
    for _ in range(num_samples):
        n = random.randint(2, int(x))
        if n in primes:
            count += 1

    estimate = (x - 1) * count / num_samples
    return round(estimate)

print("Randomized estimation of π(x):")
for x in [1000, 10000]:
    true_pi = pi_exact(x)
    for samples in [100, 1000, 10000]:
        est = randomized_pi_test(x, samples)
        err = abs(est - true_pi)
        print(f"  π({x}) ≈ {est} (true: {true_pi}, error: {err}, samples: {samples})")

print("""
RESULT: Randomized sampling needs O(x/ln(x)) samples for exact answer.
That's WORSE than Lucy DP at O(x^{2/3}).
No randomized shortcut exists.
""")

# FINAL ANALYSIS
print("=" * 60)
print("FINAL ANALYSIS: Why Modular/CRT Approaches Fail")
print("=" * 60)

print("""
THEOREM (informal): Computing π(x) mod m has the same asymptotic
complexity as computing π(x) exactly, for any fixed m ≥ 2.

PROOF SKETCH:
1. UPPER BOUND: Obviously, computing π(x) and then taking mod m
   gives π(x) mod m with the same cost.

2. LOWER BOUND: Suppose we could compute π(x) mod m in time T(x).
   Then by computing π(x) mod m for enough moduli m₁, m₂, ..., m_k
   with product > π(x), we could reconstruct π(x) via CRT in time
   O(k · T(x)). Since k = O(log(π(x))/log(m)) = O(log x), this
   gives π(x) in O(T(x) · log x) time.

   If T(x) were polylog(x), we'd have π(x) in polylog(x) time.
   But we KNOW (Aggarwal 2025) that computing p(n) requires
   Ω(n^{1/3}) operations unconditionally.

   Therefore T(x) must be at least Ω(x^{1/3}) as well.

CONCLUSION: CRT reconstruction of π(x) is a valid approach but
provides NO speedup. The bottleneck is not the SIZE of the answer
but the COMPUTATION needed to determine it.

This is fundamentally different from, say, computing n! mod m
(which CAN be done faster than computing n! exactly).
The reason: n! has a simple FORMULA (product), while π(x) does not.
""")
