#!/usr/bin/env python3
"""
PROPOSAL 17: p-adic Lifting of Prime Counting Function
═══════════════════════════════════════════════════════

IDEA: Compute pi(x) mod p for many small primes p in O(polylog) time each,
then reconstruct pi(x) exactly via CRT. The key insight is that
pi(x) mod p might be computable via:

  1. The Legendre formula: pi(x) = phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1
     Working mod p, the floor functions become modular arithmetic.

  2. For a prime q, pi(x) mod q relates to counting solutions of
     certain equations over F_q via Weil's theorem.

  3. p-adic analytic continuation: The Kubota-Leopoldt p-adic L-function
     L_p(s, chi) encodes prime-counting information mod p. If we can
     evaluate it at the right point, we get pi(x) mod p.

CONJECTURE: pi(x) mod p can be computed in O(polylog(x) * poly(p)) time
for each small prime p, and O(log(x)) primes suffice for CRT reconstruction.

If true, total cost = O(log(x)) primes * O(polylog(x)) each = O(polylog(x)).

TEST: Verify that pi(x) mod p can be computed faster than pi(x) itself
for small primes p, using the Legendre recursion mod p.
"""

import math
from functools import lru_cache

def true_nth_prime(n):
    if n < 1:
        return 2
    limit = max(100, int(n * (math.log(n) + math.log(math.log(n + 2)) + 3)))
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    primes = [i for i in range(2, limit + 1) if sieve[i]]
    if n <= len(primes):
        return primes[n - 1]
    return None

def pi_exact(x):
    """Exact prime counting function via sieve (for verification)."""
    if x < 2:
        return 0
    x = int(x)
    sieve = [True] * (x + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(x**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, x + 1, i):
                sieve[j] = False
    return sum(sieve)

def small_primes_up_to(limit):
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(2, limit + 1) if sieve[i]]

# ─── Approach A: Legendre phi mod p ───────────────────────────

def legendre_phi_mod_p(x, a, base_primes, mod_p):
    """
    Compute phi(x, a) mod p where phi is the Legendre sieve function.
    phi(x, 0) = floor(x).
    phi(x, a) = phi(x, a-1) - phi(floor(x/p_a), a-1).
    """
    if a == 0:
        return int(x) % mod_p
    if x < 1:
        return 0

    # Memoization key
    key = (int(x), a)
    if key in _phi_cache:
        return _phi_cache[key]

    val = (legendre_phi_mod_p(x, a - 1, base_primes, mod_p) -
           legendre_phi_mod_p(x // base_primes[a - 1], a - 1, base_primes, mod_p)) % mod_p

    _phi_cache[key] = val
    return val

_phi_cache = {}

def pi_mod_p_via_legendre(x, mod_p):
    """
    Compute pi(x) mod p using the Legendre/Meissel formula.
    pi(x) = phi(x, a) + a - 1 - P2(x, a) where a = pi(x^{1/3}).

    This still requires knowing base primes, so it's not truly O(polylog).
    But the modular arithmetic might allow shortcuts.
    """
    global _phi_cache
    _phi_cache = {}

    if x < 2:
        return 0

    cbrt_x = int(x ** (1.0/3.0))
    sqrt_x = int(x ** 0.5)

    # We need base primes up to cbrt(x)
    base_primes = small_primes_up_to(cbrt_x + 1)
    a = len(base_primes)

    # phi(x, a) mod p
    phi_val = legendre_phi_mod_p(x, a, base_primes, mod_p)

    # P2(x, a): sum over primes p in (cbrt_x, sqrt_x] of pi(x/p) - pi(p) + 1
    # This is the expensive part -- still O(x^{2/3}) even mod p
    mid_primes = small_primes_up_to(sqrt_x + 1)
    P2 = 0
    for p in mid_primes:
        if p <= cbrt_x:
            continue
        if p > sqrt_x:
            break
        # pi(x/p) computed mod p -- still expensive
        P2 = (P2 + pi_exact(x // p) - pi_exact(p) + 1) % mod_p

    result = (phi_val + a - 1 - P2) % mod_p
    return result

# ─── Approach B: CRT reconstruction ──────────────────────────

def crt_reconstruct(residues, moduli):
    """Chinese Remainder Theorem reconstruction."""
    M = 1
    for m in moduli:
        M *= m

    result = 0
    for ri, mi in zip(residues, moduli):
        Mi = M // mi
        # Find Mi_inv mod mi
        Mi_inv = pow(Mi, -1, mi)
        result = (result + ri * Mi * Mi_inv) % M

    return result

# ─── Approach C: Fermat quotient connection ───────────────────

def pi_mod_p_fermat(x, p):
    """
    Speculative: Can we relate pi(x) mod p to Fermat quotients?

    Wilson's theorem: (p-1)! ≡ -1 (mod p)
    Fermat quotient: q_p(a) = (a^{p-1} - 1) / p mod p

    The number of primes up to x mod p might relate to:
    sum_{prime q <= x} q_p(q) mod p

    This is speculative but testable.
    """
    primes = small_primes_up_to(int(x))
    count = len(primes) % p

    # Also compute Fermat quotient sum for comparison
    fq_sum = 0
    for q in primes:
        if q != p:
            fq = (pow(q, p - 1, p * p) - 1) // p % p
            fq_sum = (fq_sum + fq) % p

    return count, fq_sum

def run_test():
    print("PROPOSAL 17: p-adic Lifting of Prime Counting Function")
    print("=" * 60)

    # Test 1: Verify Legendre mod p gives correct residues
    print("\n--- Test 1: Legendre formula mod p correctness ---")
    test_moduli = [3, 5, 7, 11, 13, 17, 19, 23]
    test_x_values = [100, 500, 1000, 5000, 10000]

    correct = 0
    total = 0
    for x in test_x_values:
        true_pi = pi_exact(x)
        residues = []
        for p in test_moduli:
            computed = pi_mod_p_via_legendre(x, p)
            expected = true_pi % p
            match = "✓" if computed == expected else "✗"
            if computed == expected:
                correct += 1
            total += 1
            residues.append(computed)

        # CRT reconstruction
        crt_val = crt_reconstruct(residues, test_moduli)
        product = 1
        for m in test_moduli:
            product *= m
        crt_correct = (crt_val == true_pi) if true_pi < product else (crt_val == true_pi % product)

        print(f"x={x:>6}: pi(x)={true_pi:>5}, CRT={crt_val:>8}, "
              f"product={product:>10}, CRT correct: {crt_correct}")

    print(f"\nMod-p correctness: {correct}/{total}")

    # Test 2: How many primes needed for CRT?
    print("\n--- Test 2: CRT moduli needed for exact reconstruction ---")
    all_moduli = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    for x in [100, 1000, 10000, 50000, 100000]:
        true_pi = pi_exact(x)
        product = 1
        needed = 0
        for p in all_moduli:
            product *= p
            needed += 1
            if product > true_pi:
                break
        print(f"x={x:>7}: pi(x)={true_pi:>6}, need {needed} primes "
              f"(product={product:>12} > {true_pi}), "
              f"log(pi(x))={math.log(true_pi):.1f}, "
              f"ratio={needed/math.log(true_pi + 1):.2f}")

    # Test 3: Fermat quotient correlation
    print("\n--- Test 3: Fermat quotient correlation ---")
    print(f"{'p':>4} | {'x':>6} | {'pi(x) mod p':>11} | {'FQ sum mod p':>12} | {'match':>5}")
    print("-" * 55)

    matches = 0
    ftotal = 0
    for p in [3, 5, 7, 11, 13]:
        for x in [100, 500, 1000]:
            count_mod_p, fq_sum = pi_mod_p_fermat(x, p)
            match = count_mod_p == fq_sum
            if match:
                matches += 1
            ftotal += 1
            print(f"{p:>4} | {x:>6} | {count_mod_p:>11} | {fq_sum:>12} | {'✓' if match else '✗':>5}")

    print(f"\nFermat quotient matches: {matches}/{ftotal}")
    print("(Random would give ~1/p match rate)")

    # Key analysis
    print("\n--- COMPLEXITY ANALYSIS ---")
    print("Legendre mod p: Still O(x^{2/3}) per prime -- the floor divisions")
    print("  don't simplify under modular reduction because we need to know")
    print("  which integers are sieved out, not just a count mod p.")
    print()
    print("CRT approach: Needs O(log(pi(x))) = O(log(x)/log(log(x))) primes.")
    print("  If each residue took O(polylog) time, total would be O(polylog).")
    print("  BOTTLENECK: Computing pi(x) mod p still requires Legendre recursion.")
    print()
    print("VERDICT: The modular reduction does NOT shortcut the Legendre recursion.")
    print("  floor(x/p_a) mod q != floor((x mod q) / (p_a mod q)) in general.")
    print("  This is the fundamental obstacle -- floor division is not a homomorphism.")

    return correct, total

if __name__ == "__main__":
    run_test()
