"""
PROPOSAL 2: CRT Modular Reconstruction of delta(n)

IDEA: p(n) = R^{-1}(n) + delta(n), where delta(n) has O(log n) bits.
Instead of computing delta(n) directly, compute delta(n) mod q for
several small primes q, then use CRT to reconstruct delta(n).

KEY INSIGHT: delta(n) mod q = p(n) mod q - R^{-1}(n) mod q.
p(n) mod q depends on which residue class p(n) falls into mod q.

For SMALL q, there's a "prime race" among residue classes.
Chebyshev's bias tells us primes tend to favor non-residues mod q.

NOVEL TWIST: For q=2: p(n) mod 2 = 1 for n >= 2 (trivial!)
For q=3: p(n) mod 3 is 1 or 2 for n >= 3
For q=6: p(n) mod 6 is 1 or 5 for n >= 3

Can we determine p(n) mod q using the explicit formula for 
pi(x; q, a) = #{p <= x : p ≡ a mod q} with FEWER zeros?

The idea: pi(x;q,a) = li(x)/phi(q) - sum over zeros of L(s,chi) ...
Since L-functions have independent zeros, maybe the TOTAL information
needed across all residue classes is less than for pi(x) itself?
"""

import numpy as np
from sympy import prime, primepi, isprime, factorint
from sympy.ntheory import totient
from functools import reduce
import time

def primes_in_class(limit, q, a):
    """Count primes p <= limit with p ≡ a (mod q)"""
    count = 0
    for p in range(2, limit + 1):
        if isprime(p) and p % q == a:
            count += 1
    return count

def p_n_mod_q(n, q):
    """Compute p(n) mod q by finding which residue class the nth prime falls in."""
    # For small n, just compute directly
    p = prime(n)
    return p % q

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = extended_gcd(b % a, a)
    return gcd, y1 - (b // a) * x1, x1

def crt(remainders, moduli):
    """Chinese Remainder Theorem"""
    if not remainders:
        return 0, 1
    r, m = remainders[0], moduli[0]
    for i in range(1, len(remainders)):
        r2, m2 = remainders[i], moduli[i]
        g, x, _ = extended_gcd(m, m2)
        if (r2 - r) % g != 0:
            return None, None  # No solution
        lcm = m * m2 // g
        r = (r + m * ((r2 - r) // g * x % (m2 // g))) % lcm
        m = lcm
    return r, m

def crt_reconstruction(n, small_primes=[2, 3, 5, 7, 11, 13]):
    """
    Reconstruct p(n) using CRT:
    1. Compute R^{-1}(n) for approximate value
    2. Compute p(n) mod q for several small q
    3. Use CRT to narrow down delta(n)
    4. Search remaining candidates
    """
    from mpmath import mpf, log, li, mp
    mp.dps = 30
    
    t0 = time.time()
    
    # Step 1: R^{-1}(n) approximate
    x = mpf(n) * log(mpf(n))
    for _ in range(30):
        rx = li(x)
        dx = 1 / log(x)
        x = x + (mpf(n) - rx) / dx
    x_approx = int(float(x))
    
    # Step 2: Compute p(n) mod q for each q
    # (In a real algorithm, this would use L-function methods)
    # For now, compute directly to test the CRT framework
    p_actual = prime(n)
    
    remainders = []
    moduli = []
    
    for q in small_primes:
        r = p_actual % q  # In practice, would compute via L-functions
        remainders.append(r)
        moduli.append(q)
    
    # Step 3: CRT gives p(n) mod M where M = product of small_primes
    p_mod_M, M = crt(remainders, moduli)
    
    # Step 4: delta = p(n) - x_approx
    # We know delta ≡ (p_mod_M - x_approx) mod M
    delta_approx = p_actual - x_approx
    
    # How many candidates remain?
    # delta is in range [-D, D] where D ≈ sqrt(x) * log(x)
    D = int(float(x) ** 0.5 * float(log(x)))
    n_candidates = (2 * D) // M + 1
    
    # Search candidates: x_approx + p_mod_M + k*M for integer k
    candidates = []
    base = p_mod_M
    k_start = (x_approx - D - base) // M
    k_end = (x_approx + D - base) // M + 1
    for k in range(k_start, k_end + 1):
        c = base + k * M
        if c > 1 and abs(c - x_approx) <= D:
            candidates.append(c)
    
    # Filter: only primes
    prime_candidates = [c for c in candidates if isprime(c)]
    
    t1 = time.time()
    
    return {
        'x_approx': x_approx,
        'delta': delta_approx,
        'M': M,
        'p_mod_M': p_mod_M,
        'search_range': 2*D,
        'n_candidates': len(candidates),
        'n_prime_candidates': len(prime_candidates),
        'reduction_factor': (2*D) / max(1, len(candidates)),
        'time': t1 - t0,
        'correct': p_actual in prime_candidates
    }

print("=" * 70)
print("PROPOSAL 2: CRT Modular Reconstruction")
print("=" * 70)
print()

# Test with increasing number of CRT moduli
for n in [100, 1000, 5000, 10000]:
    p = prime(n)
    print(f"n={n}, p(n)={p}:")
    
    for primes_list in [[2,3], [2,3,5], [2,3,5,7], [2,3,5,7,11], [2,3,5,7,11,13]]:
        result = crt_reconstruction(n, primes_list)
        M = result['M']
        print(f"  moduli={primes_list} -> M={M:>6}, candidates={result['n_candidates']:>6}, "
              f"prime_cands={result['n_prime_candidates']:>5}, "
              f"reduction={result['reduction_factor']:.1f}x, correct={result['correct']}")
    print()

# Key question: how does the number of prime candidates scale?
print("\nScaling of prime candidates vs n (using moduli [2,3,5,7,11,13]):")
for n in [100, 200, 500, 1000, 2000, 5000, 10000]:
    result = crt_reconstruction(n, [2,3,5,7,11,13])
    print(f"  n={n:>6}: prime_candidates={result['n_prime_candidates']:>5}, "
          f"search_range={result['search_range']:>10}")

