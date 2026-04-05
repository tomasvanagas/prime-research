#!/usr/bin/env python3
"""
Multiplicative Convolution Shortcut

Key identity: The prime indicator function chi_P can be expressed via Möbius inversion:
  chi_P(n) = sum_{d|n} mu(d) * (n/d > 1 indicator)   [not quite, but related]

More precisely, the von Mangoldt function:
  Lambda(n) = -sum_{d|n} mu(d) * ln(d)

And pi(x) = sum_{n<=x} chi_P(n) where chi_P is the prime indicator.

The Dirichlet series of chi_P is:
  sum_{p prime} p^{-s} = sum_{k=1}^{inf} mu(k)/k * ln(zeta(ks))   [Euler product relation]

Idea: In Dirichlet convolution space, multiplicative functions "factorize."
If pi(x) can be expressed as a Dirichlet convolution of simpler functions
that are each O(polylog) to evaluate, the convolution itself might be fast.

Key insight from number theory: the PRIME INDICATOR is NOT multiplicative!
It's the LOG of a multiplicative function (via von Mangoldt).
This logarithmic relationship is what makes primes hard.

But what if we work with Lambda(n) instead and use partial summation?

Experiment: Test if psi(x) = sum_{n<=x} Lambda(n) has a fast Dirichlet convolution
representation that can be computed in sublinear time.

psi(x) = x - sum_rho x^rho/rho - ln(2*pi) - 0.5*ln(1-1/x^2)

The sum over zeros is again the bottleneck. But in multiplicative Fourier space,
does it become simpler?
"""

import numpy as np
from sympy import primerange, mobius, factorint, isprime, primepi, totient
from collections import defaultdict
import time
import math

def von_mangoldt(n):
    """Lambda(n) = ln(p) if n = p^k, else 0."""
    if n <= 1:
        return 0
    factors = factorint(n)
    if len(factors) == 1:
        p = list(factors.keys())[0]
        return math.log(p)
    return 0

def psi_function(x):
    """Chebyshev psi(x) = sum_{n<=x} Lambda(n)."""
    total = 0
    for n in range(2, x + 1):
        total += von_mangoldt(n)
    return total

def dirichlet_convolution(f, g, x):
    """Compute (f * g)(n) = sum_{d|n} f(d) * g(n/d) for all n <= x."""
    result = np.zeros(x + 1)
    for d in range(1, x + 1):
        for k in range(1, x // d + 1):
            result[d * k] += f(d) * g(k)
    return result

def test_mobius_decomposition(x):
    """
    Test: can we decompose the prime indicator chi_P via Mobius inversion
    into components that are each fast to sum?

    chi_P(n) = sum_{k=1}^{log_2(n)} mu(k) * chi_{P^k}(n) ... no, this is wrong.

    Actually: Lambda(n) = -sum_{d|n} mu(d) * ln(d)
    This is a Dirichlet convolution: Lambda = -mu * ln = mu * (-ln)

    So psi(x) = sum_{n<=x} Lambda(n) = sum_{n<=x} sum_{d|n} mu(d) * (-ln(d))
             = sum_{d<=x} mu(d) * (-ln(d)) * floor(x/d)

    This is a HYPERBOLA sum. The key question: can this be computed in o(x) time?

    Hyperbola method: split at sqrt(x).
    sum_{d<=sqrt(x)} mu(d)*(-ln(d)) * floor(x/d)
    + sum_{k<=sqrt(x)} (sum_{d: x/d=k} mu(d)*(-ln(d)))
    - correction

    This is essentially how Meissel-Lehmer works. Can we do better?
    """
    # Direct computation
    t0 = time.time()
    direct_sum = 0
    for d in range(1, x + 1):
        mu_d = mobius(d)
        if mu_d != 0:
            direct_sum += mu_d * (-math.log(d)) * (x // d)
    t1 = time.time()

    psi_x = psi_function(x)

    return {
        'x': x,
        'mobius_sum': direct_sum,
        'psi_x': psi_x,
        'match': abs(direct_sum - psi_x) < 0.5,
        'time': t1 - t0,
    }

def test_multiplicative_fourier(x):
    """
    Multiplicative Fourier analysis: decompose functions over (Z/nZ)*.

    For a Dirichlet character chi mod q, the "Fourier coefficient" is:
    hat{f}(chi) = sum_{n=1}^{q} f(n) * conj(chi(n))

    If the prime indicator is "sparse" in character space, we might
    compute pi(x; q, a) efficiently for many (q, a) pairs and reconstruct.

    But this is just the Dirichlet L-function approach in disguise.
    """
    # For each modulus q, decompose the prime indicator into characters
    for q in [3, 5, 7, 11]:
        primes_mod_q = defaultdict(int)
        for p in primerange(2, x + 1):
            primes_mod_q[p % q] += 1

        phi_q = int(totient(q))
        expected = primepi(x) / phi_q

        print(f"  mod {q}: distribution = {dict(primes_mod_q)}")
        print(f"    expected per class: {expected:.1f}")
        max_dev = max(abs(v - expected) for k, v in primes_mod_q.items()
                     if math.gcd(k, q) == 1)
        print(f"    max deviation: {max_dev:.1f} "
              f"({max_dev/math.sqrt(primepi(x))*math.sqrt(phi_q):.2f} sigma)")

def test_factored_pi_computation(x):
    """
    Idea: pi(x) = sum_{n<=x} chi_P(n)

    Using Mobius: chi_P(n) = sum_{d|n, d>1} mu(n/d) * (omega(d) == 1)
    ... this doesn't simplify nicely.

    Alternative: Use the identity
    ln(zeta(s)) = sum_p sum_k p^{-ks}/k = sum_p p^{-s} + O(1)

    So P(s) = sum_p p^{-s} ≈ ln(zeta(s)) - sum_p sum_{k>=2} p^{-ks}/k

    The second sum converges to O(1) for Re(s) > 1/2.

    Perron's formula: pi(x) = (1/2pi*i) integral P(s) * x^s/s ds

    This requires evaluating P(s) on a contour, which requires knowing all primes...
    CIRCULAR!

    Unless: P(s) = ln(zeta(s)) - R(s) where R(s) is the prime power correction.
    R(s) is small and convergent. So we need zeta(s) on a contour.
    This is exactly the Lagarias-Odlyzko approach.
    """
    # Test: how many terms of the prime power correction R(s) are needed?
    # R(s) = sum_{p, k>=2} 1/(k * p^{ks})
    # For s = 1, R(1) = sum_{p, k>=2} 1/(k * p^k)
    # This converges FAST because p^k grows rapidly.

    primes = list(primerange(2, 100))
    for s_val in [1.0, 0.75, 0.6, 0.51]:
        R = 0
        for p in primes:
            for k in range(2, 50):
                term = 1.0 / (k * p**(k * s_val))
                if term < 1e-15:
                    break
                R += term

        print(f"  R(s={s_val:.2f}) = {R:.6f} (converges in O(1) terms)")

    print(f"  → Prime power correction R(s) is O(1) to compute for any s > 1/2")
    print(f"  → Bottleneck is evaluating ln(zeta(s)) on a contour, which needs zeros")

def main():
    print("=" * 70)
    print("MULTIPLICATIVE CONVOLUTION SHORTCUT")
    print("=" * 70)

    # Test 1: Mobius decomposition of psi(x)
    print("\n--- Test 1: Mobius Decomposition of psi(x) ---")
    for x in [100, 500, 1000]:
        result = test_mobius_decomposition(x)
        print(f"  x={x}: mobius_sum={result['mobius_sum']:.2f}, "
              f"psi={result['psi_x']:.2f}, "
              f"match={result['match']}, time={result['time']:.3f}s")

    # Test 2: Multiplicative Fourier analysis
    print("\n--- Test 2: Multiplicative Fourier Analysis ---")
    test_multiplicative_fourier(10000)

    # Test 3: Prime power correction
    print("\n--- Test 3: Prime Power Correction Convergence ---")
    test_factored_pi_computation(10000)

    # Test 4: Hyperbola method timing
    print("\n--- Test 4: Hyperbola Method Speedup ---")
    for x in [1000, 5000, 10000]:
        # Direct: O(x)
        t0 = time.time()
        direct = psi_function(x)
        t_direct = time.time() - t0

        # Hyperbola: O(sqrt(x))
        t0 = time.time()
        sqrt_x = int(math.isqrt(x))
        # sum_{d<=sqrt(x)} mu(d)*(-ln(d)) * floor(x/d)
        part1 = 0
        for d in range(1, sqrt_x + 1):
            mu_d = mobius(d)
            if mu_d != 0:
                part1 += mu_d * (-math.log(d)) * (x // d)

        # This is only the first part; full hyperbola needs more work
        t_hyper = time.time() - t0

        print(f"  x={x}: direct={t_direct:.3f}s, partial_hyper={t_hyper:.3f}s, "
              f"speedup={t_direct/t_hyper:.1f}x")

    print("\n--- ANALYSIS ---")
    print("The multiplicative convolution approach reveals that:")
    print("1. psi(x) = -sum_{d<=x} mu(d)*ln(d)*floor(x/d) (exact)")
    print("2. This is a hyperbola sum, computable in O(sqrt(x)) time")
    print("3. But O(sqrt(x)) >> polylog for x = 10^100")
    print("4. The bottleneck is the FLOOR function: floor(x/d) introduces")
    print("   rounding that couples all prime divisors non-multiplicatively")
    print("5. In Dirichlet series space, everything factorizes beautifully")
    print("   (Euler product), but evaluating at a point x requires a")
    print("   contour integral that needs zeta zero information")
    print()
    print("VERDICT: Multiplicative structure of Dirichlet series helps")
    print("reduce to zeta evaluation, but cannot eliminate the zero sum.")
    print("The FLOOR function is the fundamental obstacle to factorization.")

if __name__ == "__main__":
    main()
