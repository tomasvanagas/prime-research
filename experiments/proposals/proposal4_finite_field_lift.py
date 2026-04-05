"""
PROPOSAL 4: Finite Field Lifting / Function Field Analogy

IDEA: Over F_q[x], the number of irreducible polynomials of degree n is:
  I_q(n) = (1/n) * sum_{d|n} mu(n/d) * q^d

This is an EXACT, CLOSED-FORM formula with O(d(n)) divisor operations.
There's no analogue of the "sqrt(x) barrier" in function fields because
the Riemann Hypothesis is PROVED (Weil, 1948).

APPROACH: Construct a sequence of function fields F_q[x] where:
1. The irreducible polynomial counting function "approximates" pi(x) in Z
2. Take a limit as q -> 1 (the "field with one element" F_1 idea)
3. Try to extract information about p(n) from the function field side

MATHEMATICAL BASIS:
- Deninger's program: Spec(Z) should be a "3-manifold" with F_1 structure
- Connes-Consani: F_1-geometry and its zeta function
- Borger: Lambda-rings as "descent data to F_1"

The F_1 zeta function of Spec(Z) should be the Riemann zeta function.
If we can compute the "irreducible counting function" over F_1 analogously
to I_q(n), we'd get pi(x) exactly.

TEST: Compare I_q(n) behavior as q->1 with actual pi(x).
"""

import numpy as np
from sympy import prime, primepi, isprime, mobius, divisors
from mpmath import mp, mpf, log, li, zeta as mpzeta
import time

mp.dps = 30

def irreducible_count_Fq(q, n):
    """Number of irreducible polynomials of degree exactly n over F_q"""
    total = 0
    for d in divisors(n):
        total += mobius(n // d) * q**d
    return total // n

def prime_count_Fq(q, n):
    """Total irreducible polynomials of degree <= n over F_q (analogue of pi(q^n))"""
    total = 0
    for k in range(1, n + 1):
        total += irreducible_count_Fq(q, k)
    return total

def F1_limit_irreducible(n):
    """
    Attempt: take q -> 1 limit of I_q(n) / q^n
    
    I_q(n) = (1/n) sum_{d|n} mu(n/d) q^d
    
    As q -> 1: I_q(n) -> (1/n) sum_{d|n} mu(n/d) = [n=1] / n
    This is 1 if n=1, 0 otherwise. Not useful.
    
    Better: consider I_q(n) / (q-1) as q -> 1 (regularized limit)
    """
    # Use L'Hopital: d/dq [sum mu(n/d) q^d] at q=1
    # = sum mu(n/d) * d * q^{d-1} |_{q=1}
    # = sum mu(n/d) * d
    # This is the Jordan totient J_1(n) = phi(n) for the sum of d*mu(n/d)
    # Wait: sum_{d|n} mu(n/d) * d = phi(n) (Euler's totient)
    
    # So regularized: I_{q->1}(n) ~ phi(n) / n * (q-1) as q -> 1
    # Not directly useful for prime counting...
    pass

def adaptive_F1_approach(x):
    """
    Alternative: for small q (q=2,3,5,...), compute pi_q(n) where q^n ≈ x,
    so n ≈ log(x)/log(q). Then extrapolate to q->1.
    
    pi_q(n) counts irreducibles of degree ≤ n over F_q.
    As q->1, n -> infinity, with q^n = x fixed.
    """
    results = {}
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        n = int(np.log(x) / np.log(q)) + 1
        pi_q = prime_count_Fq(q, n)
        # Normalized: how does pi_q(n) / (q^n / (n*log(q))) compare to pi(x)/(x/log(x))?
        if n > 0 and q**n > 1:
            normalized = pi_q / (q**n / (n * np.log(q)))
        else:
            normalized = 0
        results[q] = {'n': n, 'pi_q': pi_q, 'q^n': q**n, 'normalized': normalized}
    return results

print("=" * 70)
print("PROPOSAL 4: Finite Field Lifting")  
print("=" * 70)
print()

# Part 1: Exact formula over F_q
print("--- Part 1: Exact counting over F_q ---")
for q in [2, 3, 5, 7]:
    print(f"\nF_{q}[x]:")
    for n in range(1, 13):
        I = irreducible_count_Fq(q, n)
        pi_q = prime_count_Fq(q, n)
        total = q**n
        ratio = pi_q / total * n * np.log(q) if total > 0 else 0
        print(f"  deg≤{n:>2}: I({n})={I:>8}, pi_q({n})={pi_q:>10}, "
              f"q^n={total:>12}, pi/PNT_ratio={ratio:.4f}")

# Part 2: F_1 limit exploration
print("\n--- Part 2: q->1 extrapolation for x=1000 ---")
x = 1000
actual_pi = primepi(x)
print(f"Actual pi({x}) = {actual_pi}")

results = adaptive_F1_approach(x)
for q, data in sorted(results.items()):
    print(f"  q={q:>2}: n={data['n']:>3}, q^n={data['q^n']:>12}, "
          f"pi_q={data['pi_q']:>8}, normalized={data['normalized']:.4f}")

# Part 3: Richardson extrapolation to q=1
print("\n--- Part 3: Richardson extrapolation q -> 1 ---")
# Try to extrapolate pi_q(log_q(x)) to q=1
for x in [100, 1000, 10000]:
    actual = primepi(x)
    qs = [2, 3, 5, 7, 11, 13]
    pi_values = []
    for q in qs:
        n = max(1, int(round(np.log(x) / np.log(q))))
        pi_q = prime_count_Fq(q, n)
        pi_values.append(pi_q)
    
    # Polynomial extrapolation: fit pi_q as function of q, evaluate at q=1
    # Use polynomial in 1/q
    inv_qs = [1.0/q for q in qs]
    coeffs = np.polyfit(inv_qs, pi_values, min(4, len(qs)-1))
    extrapolated = np.polyval(coeffs, 1.0)  # q=1 means 1/q=1
    
    # Also try polynomial in q
    coeffs2 = np.polyfit(qs, pi_values, min(4, len(qs)-1))
    extrapolated2 = np.polyval(coeffs2, 1.0)
    
    print(f"  x={x:>6}: actual pi={actual:>5}, extrap(1/q)={extrapolated:>10.1f}, "
          f"extrap(q)={extrapolated2:>10.1f}")

# Part 4: New idea - use the EXACT Fq formula structure as a template
print("\n--- Part 4: Möbius inversion template ---")
print("Over F_q: I_q(n) = (1/n) * sum_{d|n} mu(n/d) * q^d")
print("Over Z:   pi(x)  = ??? similar Möbius structure ???")
print()
print("Testing: does pi(x) ≈ (1/log x) * sum_{d|N} mu(N/d) * f(x,d) ")
print("for some function f and integer N?")

# Try: pi(x) ≈ sum_{k=1}^{K} mu(k)/k * li(x^{1/k})  [This IS Riemann's R(x)]
# The analogy: k plays role of divisor, x^{1/k} plays role of q^d
# In F_q: d divides n. In Z: k ranges over 1,2,3,...

# Can we truncate to SMALL k? R(x) converges very fast!
for x in [100, 1000, 10000, 100000]:
    actual = primepi(x)
    x_mp = mpf(x)
    for K in [1, 2, 3, 5, 10]:
        R_K = mpf(0)
        for k in range(1, K+1):
            mu_k = mobius(k)
            if mu_k != 0:
                R_K += mpf(mu_k) / k * li(x_mp ** (mpf(1)/k))
        err = abs(float(R_K) - actual)
        print(f"  x={x:>7}, K={K:>2}: R_K(x)={float(R_K):>10.2f}, pi(x)={actual:>6}, error={err:>6.2f}")
    print()

