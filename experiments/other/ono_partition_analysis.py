"""
Ono-Craig-van Ittersum Partition Characterization: Computational Analysis
=========================================================================

Paper: "Integer partitions detect the primes" (PNAS 2024, arXiv:2405.06451)
Authors: William Craig, Jan-Willem van Ittersum, Ken Ono

CORRECT DEFINITIONS:
M_a(n) = sum of products of multiplicities m_1*m_2*...*m_a over all partitions
          of n using EXACTLY a distinct part sizes: n = m_1*s_1 + ... + m_a*s_a
          with 0 < s_1 < s_2 < ... < s_a and m_i >= 1.

M_1(n) = sum_{s|n} (n/s) = sigma_1(n) = sum of divisors of n
M_2(n) = sum_{s1<s2, m1*s1+m2*s2=n} m1*m2

CRITERION (Theorem 1.1): n >= 2 is prime iff
  (n^2 - 3n + 2)*M_1(n) - 8*M_2(n) = 0

QUESTION: Does this help compute pi(x) faster than O(x^{2/3})?
"""

import math
import time

def sigma1(n):
    """Sum of divisors of n."""
    s = 0
    for d in range(1, int(math.isqrt(n)) + 1):
        if n % d == 0:
            s += d
            if d != n // d:
                s += n // d
    return s

def M2(n):
    """MacMahon M_2(n) = sum m1*m2 over all n = m1*s1 + m2*s2, s1 < s2, m_i >= 1."""
    total = 0
    # Iterate over s1 from 1 upward
    for s1 in range(1, n):
        # Iterate over s2 > s1
        for s2 in range(s1 + 1, n):
            # Need m1*s1 + m2*s2 = n, m1 >= 1, m2 >= 1
            # m2*s2 = n - m1*s1, so n - m1*s1 must be positive and divisible by s2
            # m1 ranges from 1 to (n - s2) // s1
            max_m1 = (n - s2) // s1
            if max_m1 < 1:
                break  # s2 too large
            for m1 in range(1, max_m1 + 1):
                rem = n - m1 * s1
                if rem > 0 and rem % s2 == 0:
                    m2 = rem // s2
                    if m2 >= 1:
                        total += m1 * m2
        if s1 >= n - 1:
            break
    return total

# =========================================================================
# Verify criterion
# =========================================================================
def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

print("=" * 70)
print("Verifying Ono criterion: (n^2-3n+2)*sigma_1(n) - 8*M_2(n) = 0 iff prime")
print("=" * 70)

mismatches = 0
for n in range(2, 80):
    s1 = sigma1(n)
    m2 = M2(n)
    val = (n*n - 3*n + 2) * s1 - 8 * m2
    expected = is_prime(n)
    actual = (val == 0)
    if actual != expected:
        mismatches += 1
        print(f"  MISMATCH n={n}: val={val}, prime={expected}, test={actual}")
    elif n < 30 or expected:
        print(f"  n={n:>3}: sigma1={s1:>5}, M2={m2:>8}, val={val:>10}, prime={expected}")

print(f"\nMismatches: {mismatches}")

# =========================================================================
# Computational cost analysis
# =========================================================================
print("\n" + "=" * 70)
print("Computational Cost Analysis")
print("=" * 70)

# Time M2 computation for different n
for test_n in [50, 100, 200, 500]:
    t0 = time.time()
    m2 = M2(test_n)
    t1 = time.time()
    print(f"  M2({test_n}) = {m2}, time = {t1-t0:.4f}s")

print("""
COST ANALYSIS:
- M_1(n) = sigma_1(n): O(sqrt(n)) per evaluation
- M_2(n): O(n^2) per evaluation (iterate over s1, s2, m1)
- Single Ono test: O(n^2) — MUCH WORSE than trial division (O(sqrt(n)))
- pi(x) via Ono: O(x^3) — CATASTROPHICALLY WORSE than any known method

The Ono criterion uses sigma_1(n) = sum of divisors, which REQUIRES
knowing the factorization of n (or at least its divisors). This is
CIRCULAR — it tests primality using divisor structure.

M_2(n) is even more expensive: it iterates over ALL partitions of n
into exactly 2 distinct part sizes with multiplicities.
""")

# =========================================================================
# Circuit complexity of Ono criterion
# =========================================================================
print("=" * 70)
print("Circuit Complexity Analysis")
print("=" * 70)
print("""
The Ono criterion computes:
1. sigma_1(n) = sum of divisors — requires iterating divisors of n
   - Computing sigma_1(n) is in P but NOT known to be in NC or TC^0
   - In fact, sigma_1(n) is closely related to integer factoring
   - sigma_1 ∈ #P (counting divisors is a sum over O(sqrt(n)) terms)

2. M_2(n) — requires iterating over all 2-part-distinct partitions
   - Even more complex than sigma_1

3. Polynomial arithmetic and comparison to zero — in TC^0

The BOTTLENECK is computing sigma_1(n), which requires knowing divisors.
This makes the Ono criterion HARDER than BPSW (which only needs modular
exponentiation, in TC^0).

IMPORTANT: sigma_1(n) is NOT known to be in TC^0 or NC^1.
Computing sigma_1 requires factoring or equivalent.

CONCLUSION: The Ono criterion has WORSE circuit complexity than BPSW.
It DOES NOT provide a path to TC^0 primality or NC prime counting.
""")

# =========================================================================
# Generating function perspective
# =========================================================================
print("=" * 70)
print("Generating Function Perspective")
print("=" * 70)
print("""
The generating function for M_a is:
  U_a(q) = sum over 0<s1<...<sa of q^(s1+...+sa) / prod (1-q^si)^2

For a=1: U_1(q) = sum_{s>=1} q^s/(1-q^s)^2 = sum_{s>=1} sum_{m>=1} m*q^{ms}
         Coefficient of q^n = sum_{s|n} (n/s) = sigma_1(n). Correct.

To count primes using this:
  pi(x) = sum_{n=2}^{x} [(n^2-3n+2)*[q^n]U_1(q) = 8*[q^n]U_2(q)]

This is NOT a "nice" sum over generating function coefficients because
it involves COMPARING two different coefficient extractions.

No known technique converts this comparison into an efficient sum.
The indicators [condition = 0] don't factor through generating functions.

This approach provides NO shortcut for pi(x).
""")

# =========================================================================
# FINAL VERDICT
# =========================================================================
print("=" * 70)
print("FINAL VERDICT")
print("=" * 70)
print("""
The Ono-Craig-van Ittersum partition characterization (PNAS 2024):

  CLOSED for pi(x) computation.

  Failure mode: CIRCULARITY + WORSE COMPLEXITY
  1. Requires sigma_1(n) = sum of divisors, which needs factoring/divisors
  2. M_2(n) computation is O(n^2) per element — worse than trial division
  3. Circuit complexity is WORSE than BPSW (sigma_1 not known in TC^0)
  4. No generating function shortcut for counting zeros
  5. Does NOT escape the "three pillars" — it IS the prime positions pillar
     (just tests each position individually using divisor structure)

  Mathematical interest: HIGH (beautiful connection between additive and
  multiplicative number theory). Computational value: NONE.

  The key point: ANY primality characterization, no matter how elegant,
  still requires per-element evaluation to count primes. The barrier is
  not in TESTING but in COUNTING — and no per-element test, however
  clever, can beat O(x) for counting (without batch structure).
""")
