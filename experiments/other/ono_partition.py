#!/usr/bin/env python3
"""
Session 4: Ono-Craig-van Ittersum Partition-Based Prime Detection (2024)

Paper: arXiv:2405.06451, PNAS September 2024
Authors: William Craig, Jan-Willem van Ittersum, Ken Ono

CRITICAL INSIGHT (from the paper's generating functions):
  U_1(q) = sum_{n>=1} sigma_1(n) q^n
  U_2(q) = (1/8) sum_{n>=1} ((-2n+1)*sigma_1(n) + sigma_3(n)) q^n

Therefore:
  M_1(n) = sigma_1(n)   (sum of divisors, NOT number of divisors!)
  M_2(n) = ((-2n+1)*sigma_1(n) + sigma_3(n)) / 8

Main Theorem:
  An integer n >= 2 is prime if and only if:
    (n^2 - 3n + 2) * M_1(n) - 8 * M_2(n) = 0

Substituting M_2:
    (n^2 - 3n + 2) * sigma_1(n) - ((-2n+1)*sigma_1(n) + sigma_3(n)) = 0
    (n^2 - 3n + 2 + 2n - 1) * sigma_1(n) - sigma_3(n) = 0
    (n^2 - n + 1) * sigma_1(n) - sigma_3(n) = 0

So the formula reduces to:
    (n^2 - n + 1) * sigma_1(n) = sigma_3(n)

For prime p: sigma_1(p) = p+1, sigma_3(p) = p^3+1
    (p^2 - p + 1)(p + 1) = p^3 + 1  ...  both sides = (p+1)(p^2-p+1) = p^3+1  YES!

This is a KNOWN identity! Primes satisfy sigma_3(n) = (n^2-n+1)*sigma_1(n).

Follow-up papers:
  - arXiv:2409.14253: MacMahonesque functions detect cubes of primes, primes in APs
  - arXiv:2412.19180: Quasi-modularity in MacMahon partition variants
  - arXiv:2501.00580: Partition-theoretic model of prime distribution
  - arXiv:2511.20829: Computational model, exact pi(n) estimates up to 100,000
"""

import time
import math
from sympy import isprime, primepi, divisor_sigma, prime, factorint
from functools import lru_cache

print("=" * 80)
print("ONO-CRAIG-VAN ITTERSUM: INTEGER PARTITIONS DETECT PRIMES")
print("=" * 80)

# =====================================================================
# PART 1: Implement using divisor sum functions
# =====================================================================

def sigma(k, n):
    """sigma_k(n) = sum of d^k for d dividing n. O(sqrt(n))."""
    return int(divisor_sigma(n, k))

def M1(n):
    """M_1(n) = sigma_1(n) = sum of divisors of n."""
    return sigma(1, n)

def M2(n):
    """M_2(n) = ((-2n+1)*sigma_1(n) + sigma_3(n)) / 8."""
    s1 = sigma(1, n)
    s3 = sigma(3, n)
    numerator = (-2*n + 1) * s1 + s3
    assert numerator % 8 == 0, f"M_2({n}) not integer: numerator={numerator}"
    return numerator // 8

def ono_formula(n):
    """
    Compute (n^2 - 3n + 2)*M_1(n) - 8*M_2(n).
    This equals (n^2 - n + 1)*sigma_1(n) - sigma_3(n).
    """
    s1 = sigma(1, n)
    s3 = sigma(3, n)
    return (n*n - n + 1) * s1 - s3

def ono_is_prime(n):
    """Returns True if n passes the Ono-Craig partition primality criterion."""
    if n < 2:
        return False
    return ono_formula(n) == 0


# =====================================================================
# PART 2: Verify the prime-detecting formula
# =====================================================================

print("\n" + "=" * 80)
print("EXPERIMENT 1: Verify (n^2-n+1)*sigma_1(n) = sigma_3(n) iff n is prime")
print("=" * 80)

print(f"\n{'n':>4} {'prime?':>7} {'sig1':>8} {'sig3':>12} {'formula':>14} {'=0?':>5}")
print("-" * 60)

all_correct = True
for n in range(2, 51):
    s1 = sigma(1, n)
    s3 = sigma(3, n)
    val = (n*n - n + 1) * s1 - s3
    is_p = isprime(n)
    is_zero = (val == 0)
    ok = (is_zero == is_p)
    if not ok:
        all_correct = False
    flag = "" if ok else "  *** MISMATCH ***"
    if n <= 30 or not ok:
        print(f"{n:4d} {str(bool(is_p)):>7} {s1:8d} {s3:12d} {val:14d} "
              f"{'YES' if is_zero else 'no':>5}{flag}")

print(f"\nAll results correct for n=2..50: {all_correct}")

# Extended verification
print("\nExtended verification up to n=1000...")
mismatches = 0
for n in range(2, 1001):
    result = ono_is_prime(n)
    expected = bool(isprime(n))
    if result != expected:
        mismatches += 1
        print(f"  MISMATCH at n={n}: ono={result}, isprime={expected}")
if mismatches == 0:
    print("  All 999 values correct!")


# =====================================================================
# PART 3: Why the formula works -- algebraic proof
# =====================================================================

print("\n" + "=" * 80)
print("EXPERIMENT 2: Algebraic analysis of (n^2-n+1)*sigma_1(n) = sigma_3(n)")
print("=" * 80)

print("""
PROOF THAT THE FORMULA WORKS:

For n with prime factorization n = p1^a1 * p2^a2 * ... * pk^ak:
  sigma_k(n) = prod_i (p_i^{k(a_i+1)} - 1) / (p_i^k - 1)

For prime p:
  sigma_1(p) = p + 1
  sigma_3(p) = p^3 + 1 = (p+1)(p^2-p+1)
  (p^2-p+1) * sigma_1(p) = (p^2-p+1)(p+1) = p^3+1 = sigma_3(p)  CHECK!

For prime power p^a (a >= 2):
  sigma_1(p^a) = (p^{a+1}-1)/(p-1)
  sigma_3(p^a) = (p^{3(a+1)}-1)/(p^3-1)

  (p^{2a}-p^a+1) * (p^{a+1}-1)/(p-1) vs (p^{3(a+1)}-1)/(p^3-1)

  For p=2, a=2 (n=4):
    LHS = (4-2+1)*(8-1)/1 = 3*7 = 21... wait, let me recompute.
    sigma_1(4) = 1+2+4 = 7
    sigma_3(4) = 1+8+64 = 73
    (16-4+1)*7 = 13*7 = 91 != 73
    Formula value: 91 - 73 = 18 != 0.  Correct: 4 is not prime.

For n = p*q (two distinct primes):
  sigma_1(pq) = (1+p)(1+q)
  sigma_3(pq) = (1+p^3)(1+q^3)

  Formula: (p^2*q^2 - p*q + 1)(1+p)(1+q) - (1+p^3)(1+q^3)

  Expanding (1+p^3)(1+q^3) = (1+p)(1-p+p^2)(1+q)(1-q+q^2)

  So formula = (1+p)(1+q) * [(p^2*q^2-pq+1) - (1-p+p^2)(1-q+q^2)]

  Inner: p^2*q^2-pq+1 - (1-q+q^2-p+pq-p*q^2+p^2-p^2*q+p^2*q^2)
       = p^2*q^2-pq+1 - 1+q-q^2+p-pq+p*q^2-p^2+p^2*q-p^2*q^2
       = -2pq + q - q^2 + p + p*q^2 - p^2 + p^2*q
       = p(1-p) + q(1-q) + pq(p+q-2)
       = -p(p-1) - q(q-1) + pq(p+q-2)
       = p^2*q + p*q^2 - 2pq - p^2 + p - q^2 + q

  This is generically nonzero for p,q >= 2, confirming composites fail.
""")


# =====================================================================
# PART 4: Computational complexity analysis
# =====================================================================

print("=" * 80)
print("EXPERIMENT 3: Computational complexity of the Ono formula")
print("=" * 80)

print(f"\n{'n':>12} {'time(us)':>10} {'result':>8}")
print("-" * 35)

test_values = [10, 100, 1000, 10000, 100000, 1000000, 10**7, 10**8]
times = []
for n in test_values:
    t0 = time.time()
    for _ in range(10):
        result = ono_is_prime(n)
    dt = (time.time() - t0) / 10
    times.append(dt)
    print(f"{n:12d} {dt*1e6:10.1f} {str(result):>8}")

print("\nScaling analysis:")
for i in range(1, len(times)):
    if times[i-1] > 1e-7:
        ratio = math.log(times[i] / times[i-1]) / math.log(test_values[i] / test_values[i-1])
        print(f"  n={test_values[i-1]:>10d} -> {test_values[i]:>10d}: "
              f"time ratio {times[i]/times[i-1]:.1f}x, scaling ~ O(n^{ratio:.2f})")

print("""
NOTE: The dominant cost is computing sigma_1(n) and sigma_3(n), which require
factoring n. For general n, factoring is the bottleneck:
  - Trial division: O(n^{1/2})
  - Pollard rho: O(n^{1/4})
  - Number field sieve: sub-exponential in log(n)

For n ~ 10^100, factoring takes sub-exponential time -- still far too slow
for our goal, but this is AS FAST AS the Ono formula can possibly get.
The formula itself, given the factorization, is O(1).
""")


# =====================================================================
# PART 5: Can this become a COUNTING function?
# =====================================================================

print("=" * 80)
print("EXPERIMENT 4: Route from primality test to counting/finding primes")
print("=" * 80)

print("""
The Ono formula (n^2-n+1)*sigma_1(n) = sigma_3(n) is a primality criterion.

Route to pi(x):
  pi(x) = #{n <= x : (n^2-n+1)*sigma_1(n) = sigma_3(n)}

  This requires testing EACH n individually -- it is a pointwise criterion.
  Total cost: sum_{n=2}^{x} T(factor n)

  With trial division: O(x * x^{1/2}) = O(x^{3/2})
  With sieving the factorizations: O(x log x) to get all sigma values at once
  Compare: Sieve of Eratosthenes gives pi(x) in O(x log log x)

  So even using the Ono formula with pre-sieved factorizations, we gain NOTHING
  over simply sieving for primes directly.

Route to p(n) (n-th prime):
  Binary search on pi(x) = n requires evaluating pi(x) for various x.
  Best case: O(x log log x) per evaluation via sieve.
  With Meissel-Lehmer: O(x^{2/3}) without finding individual primes.

  The Ono formula offers NO advantage here -- it cannot skip the summation.

KEY INSIGHT: The formula requires FACTORING n to evaluate.
  If we could factor n, we already know if n is prime!
  The formula is circular in a computational sense:
    To test primality via Ono, we need factoring.
    But factoring already tells us if n is prime.

  The formula's value is THEORETICAL, not computational.
""")

# Demonstrate the circularity
print("Demonstrating computational circularity:")
print(f"{'n':>8} {'factor time':>14} {'ono time':>14} {'sympy isprime':>14}")
print("-" * 55)

for n in [1000003, 10000019, 100000007, 1000000007]:
    # Time factoring
    t0 = time.time()
    f = factorint(n)
    t_factor = time.time() - t0

    # Time Ono formula (includes factoring internally)
    t0 = time.time()
    r = ono_is_prime(n)
    t_ono = time.time() - t0

    # Time standard primality test
    t0 = time.time()
    r2 = bool(isprime(n))
    t_std = time.time() - t0

    print(f"{n:12d} {t_factor*1e6:12.0f}us {t_ono*1e6:12.0f}us {t_std*1e6:12.0f}us"
          f"  [both say {'prime' if r else 'composite'}]")


# =====================================================================
# PART 6: The deeper identity and what it really says
# =====================================================================

print("\n" + "=" * 80)
print("EXPERIMENT 5: The identity sigma_3(n) = (n^2-n+1)*sigma_1(n) for primes")
print("=" * 80)

print("""
The Ono-Craig-van Ittersum result, when unpacked through generating functions,
reduces to:

  n is prime <=> sigma_3(n) = (n^2 - n + 1) * sigma_1(n)

This is equivalent to saying: for the multiplicative function
  g(n) = sigma_3(n) / sigma_1(n)
we have g(n) = n^2 - n + 1 if and only if n is prime.

For prime p: g(p) = (p^3+1)/(p+1) = p^2 - p + 1.  Trivially true.

The deep content is:
  1. This extends to INFINITELY many such identities involving higher M_a
  2. The connection to partition theory and quasi-modular forms is new
  3. It shows additive structure (partitions) encodes multiplicative structure

But computationally, sigma_3(n)/sigma_1(n) = n^2-n+1 is a NUMBER THEORY
identity that requires factoring -- it cannot shortcut past factoring.
""")

# Show the ratio for various n
print(f"\n{'n':>6} {'sigma_3/sigma_1':>16} {'n^2-n+1':>10} {'match':>6} {'prime':>6}")
print("-" * 50)
for n in range(2, 31):
    s1 = sigma(1, n)
    s3 = sigma(3, n)
    ratio = s3 / s1
    expected = n*n - n + 1
    match = abs(ratio - expected) < 0.001
    print(f"{n:6d} {ratio:16.3f} {expected:10d} {'YES' if match else 'no':>6} "
          f"{'YES' if isprime(n) else '':>6}")


# =====================================================================
# PART 7: Schneider partition-theoretic pi(n) model
# =====================================================================

print("\n" + "=" * 80)
print("EXPERIMENT 6: Schneider partition-theoretic model (arXiv:2501.00580)")
print("=" * 80)

print("""
Schneider et al. propose:
  p_n ~ 1 + 2 * sum_{j=1}^{n-1} ceil(d(j)/2) + epsilon(n)

where d(j) = number of divisors.

Also from arXiv:2511.20829 (Part II): a refined computational model that is
"practically exact" up to n=100,000.

Key issues:
  1. The error term epsilon(n) is not computable without knowing primes
  2. Even the main term requires summing over all j < n
  3. This is an asymptotic MODEL, not an exact formula
  4. For p(10^100), the sum has 10^100 terms -- completely infeasible
""")

# Test the approximation
def schneider_approx_pn(n_index):
    """Approximate the n-th prime using Schneider's model (no error term)."""
    total = 1
    for j in range(1, n_index):
        d = int(divisor_sigma(j, 0))  # number of divisors
        total += 2 * math.ceil(d / 2)
    return total

print("Testing Schneider approximation p_n ~ 1 + 2*sum ceil(d(j)/2):")
print(f"{'n':>5} {'p(n)':>10} {'approx':>10} {'error':>8} {'rel%':>8}")
print("-" * 45)

for idx in [1, 2, 3, 5, 10, 25, 50, 100, 200]:
    actual = int(prime(idx))
    approx = schneider_approx_pn(idx)
    err = approx - actual
    rel = abs(err) / actual * 100
    print(f"{idx:5d} {actual:10d} {approx:10d} {err:+8d} {rel:7.2f}%")


# =====================================================================
# PART 8: Final verdict
# =====================================================================

print("\n" + "=" * 80)
print("FINAL ANALYSIS AND VERDICT")
print("=" * 80)

print("""
FINDINGS FROM THIS INVESTIGATION
=================================

1. FORMULA SIMPLIFICATION:
   The Ono-Craig-van Ittersum partition prime criterion, when the generating
   functions for M_1 and M_2 are unpacked, reduces to:

     n is prime <=> (n^2 - n + 1) * sigma_1(n) = sigma_3(n)

   This is an identity in standard divisor sum functions. The partition-theoretic
   framework provides the PROOF, but the resulting test is entirely classical.

2. COMPUTATIONAL COMPLEXITY:
   - Evaluating the formula requires computing sigma_1(n) and sigma_3(n)
   - Both require the complete factorization of n
   - Factoring is the bottleneck: O(n^{1/4}) with Pollard rho, sub-exp with NFS
   - But if we can factor n, we ALREADY know if n is prime
   - The formula is COMPUTATIONALLY CIRCULAR

3. COUNTING FUNCTION pi(x):
   - Using Ono formula: must test each n individually -> O(x * T_factor)
   - No shortcut past pointwise evaluation
   - Standard sieving is faster in every regime

4. SCHNEIDER MODEL:
   - Approximate, error term not computable
   - Requires O(n) summation even for the main term
   - No route to p(10^100)

5. FOLLOW-UP PAPERS:
   - arXiv:2409.14253: Higher MacMahonesque functions detect cube-primes and
     primes in APs, but are HARDER to compute
   - arXiv:2412.19180: Quasi-modularity explains the structure but offers no
     computational speedup
   - arXiv:2501.00580, 2511.20829: Partition model of pi(n), approximate only,
     tested up to n=100,000

VERDICT FOR p(10^100):
======================
COMPLETELY INFEASIBLE via any partition-based method.

The Ono-Craig-van Ittersum theorem is a beautiful piece of mathematics
connecting partitions to primes, but it provides ZERO computational advantage.
When reduced to explicit formulas, it becomes a classical identity in divisor
functions that requires factoring -- which already determines primality.

There is no route from partition theory to efficient prime computation.
The connection is structural/theoretical, not algorithmic.

+-------------------------------+------------------+---------------------+
| Method                        | Complexity       | p(10^100) feasible? |
+-------------------------------+------------------+---------------------+
| Ono formula (per number)      | O(n^{1/4}) factor| NO (= factoring)   |
| Ono -> pi(x) summation        | O(x^{5/4})       | NO                  |
| Schneider model               | O(n) approx      | NO                  |
| Miller-Rabin (per number)     | O(log^2 n)       | YES (test only)     |
| Meissel-Lehmer pi(x)          | O(x^{2/3})       | NO (x too large)    |
+-------------------------------+------------------+---------------------+
""")
