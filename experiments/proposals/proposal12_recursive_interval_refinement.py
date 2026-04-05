"""
PROPOSAL 4: Recursive Interval Refinement via Pratt Certificates
================================================================

IDEA: Combine three ingredients for a radically different approach:
1. R^{-1}(n) gives x_approx with |p(n) - x_approx| = O(x^{1/2}/log(x))
2. In this interval, there are O(x^{1/2}/log(x)^2) primes
3. We need to identify WHICH one is the nth

NEW INGREDIENT: Instead of sieving the entire interval, use a
HIERARCHICAL PRIMALITY CERTIFICATES approach:

Step A: Compute x_approx = R^{-1}(n) in O(polylog(n))
Step B: Determine the interval [x_approx - E, x_approx + E] where p(n) lives
Step C: Use modular arithmetic to narrow down candidates

For Step C, the key idea from ALGEBRAIC NUMBER THEORY:
- For each small prime q, we know p(n) mod q determines which residue class
  the nth prime falls in
- Primes > sqrt(x) avoid certain residue classes mod each q < sqrt(x)
- By Dirichlet's theorem, primes are ~uniform in (Z/qZ)*
- The CRT of several constraints can pinpoint p(n) in a small set

DEEPER: THE BEATTY SEQUENCE CONNECTION
- Define A = {floor(n*ln(n) + n*ln(ln(n)) + c*n) : c in [-1,2]}
- p(n) is "close to" a Beatty-like sequence with irrational slope
- The theory of Beatty sequences has efficient algorithms for the kth element

EVEN DEEPER: ITERATED LOGARITHMIC REFINEMENT
- p(n) ~ n*ln(n) + n*ln(ln(n)) - n + n*ln(ln(n))/ln(n) + ...
- Each term adds ~1 bit of precision
- After O(log n) terms, error < n^{1/2}
- Then binary search in the remaining gap

ASSUMPTION: The error in the asymptotic expansion decreases geometrically
or the modular constraints are strong enough to resolve the ambiguity.
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime, prevprime, primerange
from mpmath import mp, mpf, li, log, exp, pi as mpi, lambertw
import math

mp.dps = 50

def asymptotic_expansion_pn(n, terms=10):
    """
    Asymptotic expansion of p(n):
    p(n) = n*ln(n) + n*ln(ln(n)) - n + n*(ln(ln(n))-2)/ln(n)
           + n*(ln(ln(n))^2 - 6*ln(ln(n)) + 11)/(2*ln(n)^2) + ...

    Returns successive approximations.
    """
    if n < 6:
        return [float(prime(n))] * terms

    ln_n = math.log(n)
    lnln_n = math.log(ln_n)

    approximations = []

    # Term 0: n * ln(n)
    a0 = n * ln_n
    approximations.append(a0)

    # Term 1: + n*ln(ln(n))
    a1 = a0 + n * lnln_n
    approximations.append(a1)

    # Term 2: - n
    a2 = a1 - n
    approximations.append(a2)

    # Term 3: + n*(ln(ln(n)) - 2)/ln(n)
    a3 = a2 + n * (lnln_n - 2) / ln_n
    approximations.append(a3)

    # Term 4: + n*(ln(ln(n))^2 - 6*ln(ln(n)) + 11)/(2*ln(n)^2)
    a4 = a3 + n * (lnln_n**2 - 6*lnln_n + 11) / (2 * ln_n**2)
    approximations.append(a4)

    # Term 5: higher order
    a5 = a4 - n * (lnln_n**3 - 9*lnln_n**2 + 39*lnln_n - 53) / (6 * ln_n**3)
    approximations.append(a5)

    return approximations

def R_inv_newton(n, iterations=30):
    """R^{-1}(n) via Newton's method."""
    from sympy import mobius
    if n < 2:
        return 2.0

    x = mpf(n) * log(mpf(n))
    for _ in range(iterations):
        R_val = mpf(0)
        for k in range(1, 50):
            mu_k = mobius(k)
            if mu_k != 0:
                R_val += mpf(mu_k) / k * li(x ** (mpf(1)/k))
        R_prime = 1 / log(x)
        correction = (mpf(n) - R_val) / R_prime
        x = x + correction
        if abs(correction) < mpf('1e-15'):
            break
    return float(x)

def test_expansion_convergence():
    """Test how the asymptotic expansion converges."""
    print("Asymptotic expansion convergence:")
    print(f"{'n':>8s} {'p(n)':>10s} {'R_inv':>12s} {'err_R':>8s} | " +
          " | ".join(f"T{i}" for i in range(6)))

    for n in [100, 500, 1000, 5000, 10000]:
        p_n = int(prime(n))
        r_inv = R_inv_newton(n)
        approxs = asymptotic_expansion_pn(n, terms=6)

        errors = [abs(a - p_n) for a in approxs]
        r_err = abs(r_inv - p_n)

        err_strs = [f"{e:8.1f}" for e in errors]
        print(f"{n:8d} {p_n:10d} {r_inv:12.1f} {r_err:8.1f} | " + " | ".join(err_strs))

def interval_size_analysis():
    """
    How big is the interval we need to search after R^{-1}(n)?
    This determines the residual difficulty.
    """
    print("\nInterval size after R^{-1}(n):")
    print(f"{'n':>8s} {'p(n)':>10s} {'R_inv':>12s} {'error':>10s} {'sqrt(x)':>10s} {'ratio':>8s}")

    for n in [100, 500, 1000, 5000, 10000]:
        p_n = int(prime(n))
        r_inv = R_inv_newton(n)
        error = abs(r_inv - p_n)
        sqrtx = math.sqrt(p_n)

        print(f"{n:8d} {p_n:10d} {r_inv:12.1f} {error:10.2f} {sqrtx:10.2f} {error/sqrtx:8.4f}")

def modular_sieving_in_interval(n_target):
    """
    After locating p(n) in an interval [a, b], use modular sieving
    to eliminate candidates.

    For each small prime q, we know p(n) != 0 mod q (for q < p(n)).
    More powerfully, we can compute pi(a) mod q cheaply, which tells us
    the "offset" of p(n) within the interval modulo prime-counting mod q.
    """
    p_n = int(prime(n_target))
    r_inv = R_inv_newton(n_target)
    error = abs(r_inv - p_n)

    # Create interval
    margin = max(int(2 * error), 100)
    lo = max(2, int(r_inv - margin))
    hi = int(r_inv + margin)

    # All primes in interval
    primes_in_interval = list(primerange(lo, hi + 1))

    # How many candidates survive sieving?
    # After sieving by primes up to B, we keep floor((hi-lo)*prod((1-1/p)))
    sieving_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    density = 1.0
    for q in sieving_primes:
        density *= (1 - 1/q)

    expected_survivors = (hi - lo) * density
    actual_primes = len(primes_in_interval)

    # Find the rank of p(n) among primes in interval
    rank = primes_in_interval.index(p_n) + 1 if p_n in primes_in_interval else -1

    return {
        'n': n_target,
        'p_n': p_n,
        'interval': (lo, hi),
        'interval_size': hi - lo,
        'primes_count': actual_primes,
        'rank_in_interval': rank,
        'expected_survivors': expected_survivors,
        'prime_density': actual_primes / (hi - lo) if hi > lo else 0,
    }

def hierarchical_certificate_test(n_target):
    """
    Test a hierarchical approach:
    1. Get R^{-1}(n) -> interval of size ~sqrt(x)/log(x)
    2. Use pi(lo) to determine rank within interval
    3. Test if this can be bootstrapped recursively

    Key: pi(lo) for lo ~ n*log(n) is ANOTHER instance of the same problem.
    But smaller! pi(lo) ~ lo/log(lo) ~ n.
    So computing pi(lo) exactly requires computing p(n') for n' ~ n.
    This is CIRCULAR unless we can break the cycle.

    NON-CIRCULAR VERSION: use modular constraints.
    """
    p_n = int(prime(n_target))
    r_inv = R_inv_newton(n_target)

    # The interval
    error = p_n - r_inv
    sqrt_pn = math.sqrt(p_n)

    # If we knew pi(lo) exactly, we'd need to find the (n - pi(lo))th prime after lo
    lo = int(r_inv - sqrt_pn)
    pi_lo = int(primepi(lo))
    remaining = n_target - pi_lo  # primes to count from lo

    # How many primes do we need to enumerate?
    # We need the 'remaining'th prime after lo
    # remaining ~ n - pi(R^{-1}(n) - sqrt(x)) ~ sqrt(x)/ln(x)

    # This is still O(sqrt(x)/log(x)) primes to enumerate in the worst case.
    # But if we RECURSE: use R^{-1} to approximate pi(lo), reducing the problem.

    # Recursive depth analysis
    problem_sizes = [p_n]
    current = p_n
    for level in range(10):
        next_size = math.sqrt(current) / math.log(max(current, 3))
        problem_sizes.append(next_size)
        current = next_size
        if current < 10:
            break

    return {
        'n': n_target,
        'p_n': p_n,
        'error': error,
        'sqrt_pn': sqrt_pn,
        'remaining_primes': remaining,
        'problem_sizes': problem_sizes,
    }


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 4: Recursive Interval Refinement")
    print("=" * 80)

    test_expansion_convergence()
    interval_size_analysis()

    print("\n--- Modular sieving in interval ---")
    for n in [100, 500, 1000, 5000]:
        result = modular_sieving_in_interval(n)
        print(f"  n={n:5d}: interval=[{result['interval'][0]},{result['interval'][1]}], "
              f"size={result['interval_size']}, primes={result['primes_count']}, "
              f"rank={result['rank_in_interval']}, density={result['prime_density']:.4f}")

    print("\n--- Recursive problem size reduction ---")
    for n in [1000, 10000]:
        result = hierarchical_certificate_test(n)
        sizes = result['problem_sizes']
        print(f"  n={n}: problem sizes at each recursion level:")
        for i, s in enumerate(sizes):
            print(f"    level {i}: {s:.1f}")
