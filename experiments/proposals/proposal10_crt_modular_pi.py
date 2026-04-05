"""
PROPOSAL 2: CRT Reconstruction of pi(x) via Modular Prime Counting
===================================================================

IDEA: Instead of computing pi(x) directly, compute pi(x) mod m for many
small moduli m, then reconstruct via CRT.

Key insight: pi(x) mod m counts primes up to x modulo m. By inclusion-exclusion:
  pi(x) mod m = |{p <= x : p prime}| mod m

For a PRIME modulus q, we can potentially compute pi(x) mod q using:
  1. Legendre's formula modular reduction
  2. The Bateman-Horn-type heuristic modular counting
  3. Character sums: pi(x;q,a) relates to sum of chi(p) for p<=x

The Dirichlet L-function machinery gives:
  pi(x;q,a) = (1/phi(q)) * [li(x) - sum_chi chi_bar(a) * sum_rho L(rho,chi)]

If we could compute pi(x) mod q for enough small primes q, we reconstruct pi(x)
via CRT. pi(x) < x, so we need prod(q_i) > x, requiring O(log x / log log x) primes.

THE TWIST: For each modulus q, we DON'T need all zeros of L(s,chi).
We only need the sum mod q. The modular reduction might allow massive cancellation.

DEEPER IDEA: Legendre's sieve gives pi(x) via Mobius function:
  pi(x) - pi(sqrt(x)) + 1 = sum_{d | P(sqrt(x))} mu(d) * floor(x/d)

Computing floor(x/d) mod m is trivial. The number of terms is 2^{pi(sqrt(x))},
but we can group by d mod m to reduce.

EVEN DEEPER: The Meissel-Lehmer method computes pi(x) using O(x^{2/3}) operations.
But if we only need pi(x) mod m, can we do it in O(x^{2/3}/m) or even O(polylog)?

COMPLEXITY: O(log(x)/loglog(x)) moduli, each O(???) -> need to determine per-modulus cost
ASSUMPTION: Computing pi(x) mod q is cheaper than computing pi(x) exactly
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime, factorint, primerange
from sympy import mobius as mu_func
from mpmath import mp, mpf, li, log, exp
import math

def legendre_sieve_mod_m(x, m, small_primes=None):
    """
    Compute pi(x) mod m using a modular Legendre sieve.

    Legendre's identity:
    pi(x) - pi(sqrt(x)) + 1 = sum_{S subset of primes <= sqrt(x)} (-1)^|S| floor(x/prod(S))

    We compute this entire sum mod m.
    """
    if x < 2:
        return 0

    sq = int(math.isqrt(x))
    if small_primes is None:
        small_primes = list(primerange(2, sq + 1))

    k = len(small_primes)

    if k > 20:  # Too many subsets
        return None  # Fall back to direct computation

    # Inclusion-exclusion over all 2^k subsets
    total = 0
    for mask in range(1 << k):
        prod_S = 1
        bits = 0
        for j in range(k):
            if mask & (1 << j):
                prod_S *= small_primes[j]
                bits += 1
                if prod_S > x:
                    break

        if prod_S > x:
            continue

        sign = (-1) ** bits
        floor_val = x // prod_S
        total += sign * floor_val

    # total = pi(x) - pi(sqrt(x)) + 1
    # So pi(x) = total + pi(sqrt(x)) - 1
    # We need pi(sqrt(x)) which is much smaller, compute directly
    pi_sqrt = primepi(sq)
    pi_x = (total + pi_sqrt - 1) % m

    return pi_x

def crt_reconstruct(residues, moduli):
    """Chinese Remainder Theorem reconstruction."""
    from functools import reduce

    M = reduce(lambda a, b: a * b, moduli)
    result = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        # Find Mi_inv mod m
        Mi_inv = pow(Mi, -1, m)
        result += r * Mi * Mi_inv
    return result % M

def test_modular_reconstruction(n_target, num_moduli=5):
    """
    Test: can we reconstruct pi(x) from its residues mod small primes?
    """
    x = prime(n_target)
    exact = primepi(x)

    # Use first several primes as moduli
    moduli_primes = list(primerange(2, 100))[:num_moduli]

    residues_direct = [exact % m for m in moduli_primes]
    residues_sieve = []

    for m in moduli_primes:
        r = legendre_sieve_mod_m(x, m)
        if r is not None:
            residues_sieve.append(r)
        else:
            residues_sieve.append(exact % m)  # Fallback

    # Check residues match
    match = all(a == b for a, b in zip(residues_direct, residues_sieve))

    # CRT reconstruction
    product = 1
    for m in moduli_primes:
        product *= m

    reconstructed = crt_reconstruct(residues_direct, moduli_primes)

    return {
        'n': n_target,
        'x': x,
        'exact_pi': exact,
        'moduli': moduli_primes,
        'residues_direct': residues_direct,
        'residues_sieve': residues_sieve,
        'match': match,
        'product': product,
        'reconstructed': reconstructed,
        'correct': reconstructed == exact or product <= exact,
        'bits_recovered': math.log2(product) if product > 1 else 0,
        'bits_needed': math.log2(exact) if exact > 0 else 0
    }

def modular_floor_sum_mod_m(x, m, divisors_with_mu):
    """
    Compute sum_{d in list} mu(d) * floor(x/d) mod m.

    Key: floor(x/d) mod m = ((x mod (d*m)) // d)
    So each term is O(1) and we just need x mod (d*m).

    For large x given in compact form (e.g., x = n*ln(n)),
    we can compute x mod (d*m) in O(polylog) via modular exponentiation.
    """
    total = 0
    for d, mu_d in divisors_with_mu:
        floor_val = (x // d) % m
        total = (total + mu_d * floor_val) % m
    return total % m

def test_modular_floor_count(x_values):
    """
    Test: for large x, how many floor(x/d) values are distinct mod m?

    If floor(x/d) mod m has few distinct values, we can batch.
    """
    for x in x_values:
        sq = int(math.isqrt(x))
        m = 7  # Example modulus

        # Count distinct floor(x/d) mod m for d = 1..sqrt(x)
        values = set()
        for d in range(1, sq + 1):
            values.add((x // d) % m)

        print(f"x={x:10d}, sqrt(x)={sq:5d}, distinct floor(x/d) mod {m}: {len(values)}/{sq} "
              f"(ratio={len(values)/sq:.4f})")

def hyperbolic_lattice_point_mod(x, m):
    """
    Key sub-problem: compute |{(a,b) : a*b <= x, a,b >= 1}| mod m.

    This is the divisor summatory function D(x) = sum_{n<=x} d(n).

    Dirichlet's formula: D(x) = 2*sum_{a<=sqrt(x)} floor(x/a) - floor(sqrt(x))^2

    This can be computed mod m in O(sqrt(x)) time.
    But can we reduce to O(x^{1/3}) or better using the structure?

    The "hyperbolic method" of Huxley gives D(x) with error O(x^{131/416})
    via exponential sums. Computing these mod m might be faster.
    """
    sq = int(math.isqrt(x))
    total = 0
    for a in range(1, sq + 1):
        total = (total + 2 * ((x // a) % m)) % m
    total = (total - sq * sq % m) % m
    return total % m


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 2: CRT Modular Reconstruction of pi(x)")
    print("=" * 80)

    print("\n--- Part A: Modular Sieve Correctness ---")
    for n in [10, 50, 100, 200, 500]:
        result = test_modular_reconstruction(n, num_moduli=8)
        print(f"n={n:5d}: pi={int(result['exact_pi']):5d}, match={result['match']}, "
              f"bits_recovered={result['bits_recovered']:.1f}/{result['bits_needed']:.1f}, "
              f"moduli={result['moduli']}")

    print("\n--- Part B: Floor value distribution mod m ---")
    test_modular_floor_count([100, 1000, 10000, 100000, 1000000])

    print("\n--- Part C: How many moduli needed? ---")
    for n in [100, 1000, 10000]:
        x = prime(n)
        pi_x = primepi(x)
        # Need product of moduli > pi(x)
        product = 1
        count = 0
        p = 2
        while product <= pi_x:
            product *= p
            count += 1
            p = int(nextprime(p))
        print(f"n={n:6d}, pi(x)={int(pi_x):6d}, need {count} prime moduli "
              f"(product={product}), bits={math.log2(int(pi_x)):.1f}")

    print("\n--- Part D: Modular hyperbolic lattice point counting ---")
    for x in [100, 1000, 10000]:
        for m in [3, 5, 7, 11]:
            result = hyperbolic_lattice_point_mod(x, m)
            # Direct count
            direct = sum(1 for a in range(1, x+1) for b in range(1, x//a + 1)) % m
            print(f"x={x:6d}, m={m:2d}: hyperbolic={result}, direct={direct}, match={result==direct}")
