#!/usr/bin/env python3
"""
Commutative Matrix Powering: The AKS Special Case
===================================================
Session 11 experiment.

KEY OBSERVATION: The matrix M in our AKS reduction represents multiplication
by (x+a) in the COMMUTATIVE ring Z_n[x]/(x^r - 1).

For GENERAL matrix powering (M^n mod m), the best known is NC^1 (not TC^0).
But our matrix is a COMPANION MATRIX of a commutative ring.

Question: Does commutativity help put matrix powering in TC^0?

For SCALAR powering (b^n mod m), the Allender (1999) TC^0 algorithm works by:
1. Factor m = p1^a1 * ... * pk^ak (via CRT, factoring is TC^0 for poly-size m)
2. For each pi^ai: compute b^n mod pi^ai using:
   - If gcd(b, pi) = 1: use b^{phi(pi^ai)} = 1 mod pi^ai
   - Compute n mod phi(pi^ai), then b^{n mod phi} mod pi^ai
   - The exponent n mod phi is a DIVISION problem (TC^0 by HAB)
   - The base-to-reduced-exponent is an iterated multiplication (TC^0 by HAB)
3. CRT to combine (TC^0)

For COMMUTATIVE RING powering (element^n in Z_n[x]/(x^r-1)):
- The ring R = Z_n[x]/(x^r-1) is commutative
- An element a in R has a^{|R*|} = 1 if a is a unit
- |R*| depends on the structure of Z_n[x]/(x^r-1)
- If we knew |R*|, we could reduce n mod |R*| and then compute a^{n mod |R*|}

But |R*| depends on n (since the ring is Z_n[x]/(x^r-1))!
For prime n: Z_n is a field, and x^r - 1 factors into cyclotomic polynomials
over Z_n. The structure depends on the order of n mod r.

This makes the "reduce exponent" step circular for AKS
(we're trying to determine if n is prime!).

Let's investigate computationally.
"""

import math
from functools import reduce


def poly_mult_mod(p1, p2, r, n):
    """Multiply two polynomials in Z_n[x]/(x^r - 1)."""
    result = [0] * r
    for i in range(len(p1)):
        if p1[i] == 0:
            continue
        for j in range(len(p2)):
            if p2[j] == 0:
                continue
            idx = (i + j) % r
            result[idx] = (result[idx] + p1[i] * p2[j]) % n
    return result


def poly_pow_mod(base, exp, r, n):
    """Compute base^exp in Z_n[x]/(x^r - 1) via repeated squaring."""
    result = [0] * r
    result[0] = 1
    b = list(base)
    while exp > 0:
        if exp & 1:
            result = poly_mult_mod(result, b, r, n)
        b = poly_mult_mod(b, b, r, n)
        exp >>= 1
    return result


def multiplicative_order(a, m):
    """Compute ord_m(a) = smallest k > 0 with a^k = 1 mod m."""
    if math.gcd(a, m) != 1:
        return None
    k = 1
    current = a % m
    while current != 1:
        current = (current * a) % m
        k += 1
        if k > m:
            return None
    return k


def analyze_ring_structure(n, r):
    """
    Analyze the structure of R = Z_n[x]/(x^r - 1).

    If n is prime:
      Z_n is a field, and x^r - 1 factors into cyclotomic polynomials
      over Z_n. The factorization depends on ord_r(n).

    The multiplicative group R* has order related to n^{ord_r(n)} - 1.
    """
    print(f"\n  Ring R = Z_{n}[x]/(x^{r} - 1):")
    print(f"    Ring size: n^r = {n}^{r} = {n**r}")

    is_n_prime = all(n % i != 0 for i in range(2, min(n, 10000)))
    print(f"    n is {'prime' if is_n_prime else 'composite'}")

    if is_n_prime and r < n:
        ord_r_n = multiplicative_order(n, r)
        if ord_r_n:
            print(f"    ord_{r}({n}) = {ord_r_n}")
            # Number of irreducible factors of x^r - 1 over Z_n
            # = number of cyclotomic polynomials Phi_d(x) for d | r,
            #   each splitting into phi(d)/ord_d(n) factors of degree ord_d(n)
            # Total: r/ord_r_n irreducible factors... approximately

            # For the special case gcd(r, n) = 1:
            # x^r - 1 = product over d|r of Phi_d(x)
            # Each Phi_d(x) splits into phi(d)/f_d factors over Z_n,
            # where f_d = ord_d(n)

            # The ring decomposes as product of fields:
            # R ~ product of GF(n^{f_d})
            # for each factor

            print(f"    Approximate decomposition: R ~ product of GF({n}^{ord_r_n})")
            print(f"    |R*| approximately: ({n}^{ord_r_n} - 1)^{r // ord_r_n}")
        else:
            print(f"    gcd({n}, {r}) > 1")


def test_reduced_exponent_approach():
    """
    Test: can we compute (x+a)^n mod (x^r-1, n) using reduced exponents?

    For a prime n and the element g = (x+a) in R = Z_n[x]/(x^r-1):
    The order of g divides |R*|.
    If we knew ord(g), we could compute g^{n mod ord(g)} instead of g^n.

    The question: is ord(g) easy to compute?
    """
    print(f"\n{'='*70}")
    print("TEST: Reduced Exponent Approach for AKS")
    print(f"{'='*70}")

    test_cases = [
        (7, 5, 1),
        (11, 7, 1),
        (13, 11, 1),
        (17, 13, 1),
        (23, 19, 1),
        (29, 23, 1),
        (31, 29, 1),
    ]

    for n, r, a in test_cases:
        # Compute (x+a)^n directly
        base = [0] * r
        base[0] = a % n
        if r > 1:
            base[1] = 1

        result_direct = poly_pow_mod(base, n, r, n)

        # Find the multiplicative order of (x+a) in the ring
        current = list(base)
        identity = [0] * r
        identity[0] = 1

        order = None
        elem = list(base)
        for k in range(1, n**r):
            if elem == identity:
                order = k
                break
            elem = poly_mult_mod(elem, base, r, n)
            if k > 1000:  # safety limit
                break

        analyze_ring_structure(n, r)

        if order:
            reduced_exp = n % order
            result_reduced = poly_pow_mod(base, reduced_exp, r, n)
            match = (result_direct == result_reduced)
            print(f"    n={n}, r={r}: ord(x+{a}) = {order}, "
                  f"n mod ord = {reduced_exp}, match: {match}")

            # The key question: how big is ord(g)?
            # If ord(g) = O(polylog(n)), the reduction helps!
            # If ord(g) = Theta(n), it doesn't help.
            print(f"    ord/n ratio: {order/n:.2f}")
        else:
            print(f"    n={n}, r={r}: order > 1000 (too large to find)")


def analyze_crt_decomposition():
    """
    Analyze: Can the CRT decomposition of the ring help?

    If x^r - 1 = f_1(x) * ... * f_s(x) mod n, then:
    R = Z_n[x]/(x^r-1) ~ Z_n[x]/(f_1) x ... x Z_n[x]/(f_s)

    In each factor, (x+a) maps to some element. Computing (x+a)^n in each
    factor independently, then combining via CRT, might be more efficient
    if the factors have special structure.
    """
    print(f"\n{'='*70}")
    print("ANALYSIS: CRT Decomposition of the Polynomial Ring")
    print(f"{'='*70}")

    print("""
  For prime n and r with gcd(n,r) = 1:
    x^r - 1 = prod_{d|r} Phi_d(x)  mod n

  Each Phi_d(x) splits into phi(d)/f_d irreducible factors mod n,
  where f_d = ord_d(n) = multiplicative order of n in (Z/dZ)*.

  The ring decomposes:
    R ~ prod_{d|r} prod_{j=1}^{phi(d)/f_d} GF(n^{f_d})

  In each GF(n^{f_d}), the element (x+a) has order dividing n^{f_d} - 1.

  To compute (x+a)^n in GF(n^{f_d}):
  - If f_d = 1: we're working in Z_n, and scalar powering IS in TC^0.
  - If f_d > 1: we need powering in GF(n^{f_d}), which is a FIELD.
    Powering in a field is like scalar powering: reduce n mod (n^{f_d} - 1).

  CRUCIAL POINT: In a FINITE FIELD GF(q), the Frobenius endomorphism gives:
    a^q = a for all a in GF(q)

  So (x+a)^n in GF(n^{f_d}) depends only on n mod (n^{f_d} - 1).

  For the AKS test: we need (x+a)^n mod (x^r-1, n).
  If n IS prime: x^r-1 splits into factors of degree f_d.
    In each factor, a^n = a (Frobenius). So the test is trivially satisfied!
    This is CORRECT -- AKS uses this property: if n is prime, the test passes.

  The AKS test works by checking that THE CONVERSE holds:
  if (x+a)^n = x^n + a for ENOUGH values of a, then n must be prime.

  For COMPOSITE n: the ring Z_n is NOT a field. The decomposition doesn't
  give nice fields. So the test FAILS for composites (some a gives wrong result).

  COMPUTATIONAL IMPLICATION:
  Computing (x+a)^n in each CRT component is equivalent to computing the
  Frobenius map on the quotient ring. For prime n, this is trivial.
  For composite n, this is where the computational difficulty lies.

  But we DON'T KNOW if n is prime -- that's what we're testing!
  So we can't use the "assume n is prime and apply Frobenius" shortcut.

  This is the FUNDAMENTAL CIRCULARITY: to use the algebraic structure of
  Z_n[x]/(x^r-1) to speed up the computation, we'd need to know if n is
  prime, but that's what we're computing.
""")

    print("CONCLUSION:")
    print("  The CRT decomposition of Z_n[x]/(x^r-1) does NOT help because:")
    print("  1. For prime n: the test is trivially satisfied (Frobenius)")
    print("  2. For composite n: the ring structure is complicated")
    print("  3. We don't know which case applies (that's the question!)")
    print("  4. This is a deep form of circularity specific to primality testing")
    print()
    print("  The matrix powering M^n mod m question remains OPEN.")
    print("  Commutativity of the polynomial ring does NOT obviously help")
    print("  because the useful properties (Frobenius, field structure)")
    print("  only hold for prime n, and primality is what we're testing.")


def main():
    test_reduced_exponent_approach()
    analyze_crt_decomposition()

    print(f"\n{'='*70}")
    print("FINAL ASSESSMENT")
    print(f"{'='*70}")
    print("""
  The TC^0 status of matrix powering M^n mod m (for polylog-dimensional M)
  remains OPEN, even for the special case of companion matrices in
  commutative polynomial rings.

  The commutativity does NOT help because:
  - The "nice" algebraic properties depend on n being prime (circular!)
  - For composite n, the ring has zero divisors and the structure is opaque
  - We can't assume primality to speed up the primality test

  However, this analysis reveals that the problem has TWO subcases:
  1. n prime: computation is trivial (Frobenius theorem)
  2. n composite: computation seems hard (no field structure)

  A potential approach: show that ANY efficient algorithm for case 2
  would constitute a primality test, making the problem AT LEAST as hard
  as primality testing. This would be a CONDITIONAL lower bound.

  STATUS: OPEN. The matrix powering reduction to PRIMES in TC^0 stands.
  Commutativity explored but does not help due to circularity.
""")


if __name__ == "__main__":
    main()
