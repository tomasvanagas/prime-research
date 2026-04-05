"""
CRT Prime Locator: Fresh approach to computing p(n) via modular arithmetic.

IDEA: Instead of computing p(n) directly, compute p(n) mod m for many moduli m,
then reconstruct via CRT.

Key question: Can we compute p(n) mod q efficiently?

p(n) mod q relates to counting primes in arithmetic progressions.
If we know π(x; q, a) for all a mod q, we can determine which residue class
the nth prime falls in.

Approach:
1. Use R^{-1}(n) to get x_approx (polylog time, ~50% digits correct)
2. For each small prime q, compute π(x_approx; q, a) for each a
3. This tells us p(n) mod q
4. CRT reconstruct the exact answer

The question is whether step 2 can be done in polylog time.
"""

import sympy
from sympy import primepi, prime, nextprime
from collections import defaultdict
import time
import math

def li(x):
    """Logarithmic integral via series"""
    if x <= 1:
        return 0
    from mpmath import mp, li as mpli
    mp.dps = 50
    return float(mpli(x))

def R_inverse_approx(n):
    """Approximate nth prime via inverse Riemann R function"""
    # Simple approximation: p(n) ≈ n * ln(n) for large n
    if n < 6:
        return [2, 3, 5, 7, 11][n-1]
    x = n * math.log(n) + n * math.log(math.log(n))
    # Newton iteration on R(x) = n
    for _ in range(20):
        rx = li(x) - 0.5 * li(x**0.5)
        drx = 1.0 / math.log(x) if x > 1 else 1
        if abs(drx) < 1e-30:
            break
        x = x - (rx - n) / drx
        if x < 2:
            x = 2
    return x

def count_primes_in_progression(x, q, a):
    """Count primes ≤ x that are ≡ a (mod q). Brute force for testing."""
    count = 0
    p = 2
    while p <= x:
        if p % q == a % q:
            count += 1
        p = sympy.nextprime(p)
    return count

def crt_reconstruct(residues, moduli):
    """Chinese Remainder Theorem reconstruction"""
    from functools import reduce
    M = reduce(lambda a, b: a * b, moduli)
    result = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        # Find inverse of Mi mod m
        yi = pow(Mi, -1, m)
        result += r * Mi * yi
    return result % M

def test_crt_approach(n, moduli=[3, 5, 7, 11, 13]):
    """
    Test: can we determine p(n) mod q from prime counting in progressions?
    """
    target = sympy.prime(n)
    x_approx = R_inverse_approx(n)

    print(f"\np({n}) = {target}")
    print(f"R^{{-1}}({n}) ≈ {x_approx:.1f}")
    print(f"Approximation error: {abs(x_approx - target):.1f}")

    residues = []
    for q in moduli:
        r = target % q
        residues.append(r)
        print(f"  p({n}) mod {q} = {r}")

    # Can we determine these residues without knowing target?
    # We need: among primes near x_approx, which residue class does the nth one fall in?

    # Reconstruct
    reconstructed = crt_reconstruct(residues, moduli)
    M = 1
    for m in moduli:
        M *= m

    print(f"\nCRT product M = {M}")
    print(f"CRT reconstruction: {reconstructed}")
    print(f"Actual mod M: {target % M}")
    print(f"Match: {reconstructed == target % M}")

    # How many moduli do we need?
    # For p(n) with d digits, we need M > p(n), so we need enough moduli
    # that their product exceeds p(n).
    # For p(10^100), p(n) ≈ 10^100 * ln(10^100) ≈ 2.3 * 10^102
    # Product of first k primes: primorial(k)
    # primorial(35) ≈ 10^36... we need primorial covering 10^102
    # That's about 240 primes (by prime number theorem, sum of logs ≈ 240)

    return reconstructed == target % M

def test_progression_counting_structure(max_x=1000):
    """
    Investigate: is π(x; q, a) easier to compute than π(x)?

    By Dirichlet, π(x; q, a) ≈ Li(x)/φ(q) for (a,q)=1.
    The error term involves zeros of Dirichlet L-functions.

    KEY QUESTION: Is the error in π(x; q, a) more structured/compressible
    than the error in π(x)?
    """
    from sympy import totient

    print("\n=== Structure of π(x; q, a) errors ===")
    for q in [3, 5, 7, 11]:
        phi_q = totient(q)
        print(f"\nq = {q}, φ(q) = {phi_q}")

        # Compute errors for each residue class
        for a in range(1, q):
            if math.gcd(a, q) != 1:
                continue

            counts = []
            expected = []
            errors = []

            for x in range(10, max_x + 1, 10):
                cnt = count_primes_in_progression(x, q, a)
                exp = li(x) / phi_q
                counts.append(cnt)
                expected.append(exp)
                errors.append(cnt - exp)

            # Analyze error structure
            if errors:
                max_err = max(abs(e) for e in errors)
                avg_err = sum(abs(e) for e in errors) / len(errors)
                print(f"  a={a}: max_error={max_err:.2f}, avg_error={avg_err:.2f}")

def test_idea_crt_with_smooth_base(max_n=100):
    """
    REFINED IDEA:
    1. Compute x_approx = R^{-1}(n) (polylog, gives ~50% digits)
    2. The true p(n) is within some window of x_approx
    3. Within that window, determine p(n) mod q for small q
    4. This narrows p(n) to a smaller set of candidates
    5. CRT gives exact answer if window < product of moduli

    The key question: how wide is the window?
    By PNT, consecutive primes near x have gap ≈ ln(x).
    But R^{-1}(n) error is much larger: O(x^{1/2}).
    """
    print("\n=== CRT with smooth base: window analysis ===")

    for n in [100, 1000, 5000, 10000]:
        target = sympy.prime(n)
        x_approx = R_inverse_approx(n)
        error = abs(x_approx - target)
        rel_error = error / target

        # How many CRT moduli to cover the window?
        window = 2 * error + 100  # conservative

        # Primorial needed
        primorial = 1
        p = 2
        k = 0
        while primorial < window:
            primorial *= p
            p = sympy.nextprime(p)
            k += 1

        print(f"n={n}: p(n)={target}, approx={x_approx:.0f}, "
              f"error={error:.0f} ({rel_error:.4%}), "
              f"window={window:.0f}, need {k} CRT moduli")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: CRT Prime Locator")
    print("=" * 60)

    # Test basic CRT approach
    for n in [10, 50, 100, 500]:
        test_crt_approach(n)

    # Test window analysis
    test_idea_crt_with_smooth_base()

    # Structure analysis (small scale due to brute force)
    test_progression_counting_structure(200)
