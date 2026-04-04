#!/usr/bin/env python3
"""
Can Helfgott-Thompson's O(x^{3/5}) technique transfer to pi(x)?
================================================================
Session 11 experiment.

H-T (arXiv:2101.08773) compute M(x) = sum mu(n) in O(x^{3/5+eps}) by:
1. Decompose: M(x) = Q_even(x) - Q_odd(x) (even/odd squarefree counts)
2. Sieve by small primes up to z = x^{1/5}
3. For each small-prime pattern, count large-prime contributions
4. Large primes have at most 4 factors (since z^5 = x), so count k-almost-primes

KEY QUESTION: Can we do the same for pi(x)?

pi(x) counts numbers with EXACTLY 1 prime factor (the primes themselves).
This is a special case of "k-almost-prime counting" with k=1.

The H-T approach works for M(x) because:
- mu(n) involves ALL squarefree numbers (many factorization patterns)
- Cancellation between even and odd patterns reduces the sum
- Small prime patterns enumerate efficiently via Buchstab

For pi(x):
- We only count numbers with omega(n) = 1 (primes)
- There's NO cancellation to exploit (positive count)
- The "small prime pattern" is empty (primes have no small factors except themselves)

Let's verify this analysis computationally and explore edge cases.
"""

import math
import time
from itertools import combinations
from functools import lru_cache


def sieve_smallest_factor(limit):
    """Compute smallest prime factor for each n up to limit."""
    spf = list(range(limit + 1))
    for i in range(2, int(limit**0.5) + 1):
        if spf[i] == i:  # i is prime
            for j in range(i*i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factorize(n, spf):
    """Factorize n using smallest prime factor table."""
    factors = []
    while n > 1:
        p = spf[n]
        exp = 0
        while n % p == 0:
            n //= p
            exp += 1
        factors.append((p, exp))
    return factors


def omega(n, spf):
    """Number of distinct prime factors."""
    return len(factorize(n, spf))


def is_squarefree(n, spf):
    """Check if n is squarefree."""
    return all(e == 1 for _, e in factorize(n, spf))


def analyze_ht_decomposition(x):
    """
    Analyze how H-T decomposes M(x) and why it doesn't transfer to pi(x).
    """
    print(f"\n{'='*70}")
    print(f"H-T Decomposition Analysis for x = {x}")
    print(f"{'='*70}")

    spf = sieve_smallest_factor(x)

    # Classify all integers 1..x by their factorization properties
    primes = []
    sf_even = []  # squarefree, even omega
    sf_odd = []   # squarefree, odd omega
    non_sf = []   # not squarefree

    for n in range(2, x + 1):
        if spf[n] == n:
            primes.append(n)
        facts = factorize(n, spf)
        if all(e == 1 for _, e in facts):
            if len(facts) % 2 == 0:
                sf_even.append(n)
            else:
                sf_odd.append(n)
        else:
            non_sf.append(n)

    # Note: 1 is squarefree with 0 prime factors (even)
    M_x = len(sf_even) + 1 - len(sf_odd)  # +1 for n=1

    print(f"  Primes (pi(x)):        {len(primes)}")
    print(f"  Squarefree, even omega: {len(sf_even)} (+ 1 for n=1)")
    print(f"  Squarefree, odd omega:  {len(sf_odd)}")
    print(f"  Non-squarefree:         {len(non_sf)}")
    print(f"  M(x) = {len(sf_even)+1} - {len(sf_odd)} = {M_x}")

    # H-T decomposition: sieve by primes up to z = x^{1/5}
    z = int(x ** 0.2) + 1
    small_primes = [p for p in range(2, z + 1) if spf[p] == p]
    print(f"\n  Sieve level z = x^(1/5) = {z}")
    print(f"  Small primes: {small_primes}")

    # For M(x): count by small-prime pattern
    # A "pattern" S is a subset of small primes.
    # For each pattern S, count n <= x such that:
    #   - n is squarefree
    #   - n is divisible by exactly the primes in S (among small primes)
    #   - all other prime factors of n are > z
    #   - sign contribution: (-1)^{omega(n)}

    # For each pattern S with product P_S:
    #   contribution = (-1)^{|S|} * sum_{m <= x/P_S, gcd(m, P_small)=1, squarefree}
    #                              (-1)^{omega(m)}
    # where the sum over m counts squarefree numbers coprime to all small primes

    # These m's have all prime factors > z, so omega(m) <= 4 (since z^5 >= x)
    # For fixed omega(m) = k, counting such m is a k-almost-prime counting problem

    # Now for pi(x):
    # pi(x) = #{n <= x : omega(n) = 1, Omega(n) = 1}
    # Using the H-T sieve-by-small-primes:
    # Small primes themselves contribute pi(z) primes
    # Large primes (p > z) contribute pi(x) - pi(z)
    # NO pattern decomposition helps because primes have NO small prime factors

    print(f"\n  For pi(x):")
    print(f"  Small primes (p <= z): pi(z) = {len(small_primes)}")
    print(f"  Large primes (p > z):  {len(primes) - len(small_primes)}")
    print(f"  Total pi(x) = {len(primes)}")

    print(f"\n  H-T pattern decomposition for M(x):")
    pattern_count = 0
    for k in range(len(small_primes) + 1):
        for S in combinations(small_primes, k):
            P_S = 1
            for p in S:
                P_S *= p
            if P_S <= x:
                # Count m <= x/P_S with all factors > z, squarefree
                bound = x // P_S
                count = 0
                for m in range(1, bound + 1):
                    facts_m = factorize(m, spf) if m > 1 else []
                    if all(p > z for p, _ in facts_m) and all(e == 1 for _, e in facts_m):
                        count += 1
                if count > 0:
                    sign = (-1) ** k
                    pattern_count += 1
                    if k <= 2:  # only print first few
                        print(f"    Pattern {S}: P_S={P_S}, bound={bound}, "
                              f"count(large sf)={count}, sign={sign:+d}")
    print(f"  Total patterns with contribution: {pattern_count}")
    print(f"  Total small-prime subsets: {2**len(small_primes)}")

    # KEY ANALYSIS: For pi(x), the "pattern" is always empty (primes have
    # no small prime factors other than themselves). So:
    # - Pattern {} (empty): counts numbers <= x with all factors > z
    #   This includes 1-almost-primes (primes), 2-almost-primes (p*q),
    #   3-almost-primes (p*q*r), and 4-almost-primes (p*q*r*s)
    #   where all p,q,r,s > z.
    # - Pattern {p} for small p: counts numbers <= x/p with all factors > z
    #   These are NOT primes (they're p times something)

    # So pi(x) = pi(z) + (count of primes > z and <= x)
    #          = pi(z) + A_1(x, z)
    # where A_k(y, z) = #{n <= y : n has exactly k prime factors, all > z}

    # A_1(x, z) is counting primes in (z, x]. This IS pi(x) - pi(z).
    # Computing A_1 efficiently is THE SAME PROBLEM as computing pi(x).

    # For M(x), H-T compute A_0 + A_2 + A_4 - A_1 - A_3 (with adjustments).
    # The key: A_0 = 1 (trivial), A_2, A_3, A_4 can be expressed as sums of
    # pi-like functions at smaller arguments. And with cancellation, the total
    # work is O(x^{3/5}).

    # But for pi(x), we ONLY need A_1, which is irreducible.

    print(f"\n  KEY ANALYSIS:")
    print(f"  M(x) = sum_k (-1)^k * [sum over patterns of A_k terms]")
    print(f"       = A_0 - A_1 + A_2 - A_3 + A_4 (large factors only)")
    print(f"       + contributions from patterns with small primes")
    print(f"  ")
    print(f"  pi(x) = pi(z) + A_1(x, z)")
    print(f"  where A_1(x, z) = count of primes in (z, x]")
    print(f"  ")
    print(f"  H-T exploit cancellation: A_0 - A_1 + A_2 - A_3 + A_4")
    print(f"  For M(x), the alternating sum has cancellation that reduces work.")
    print(f"  For pi(x), we need A_1 ALONE -- no cancellation possible.")

    # Verify A_k counts
    A = [0, 0, 0, 0, 0]
    for n in range(1, x + 1):
        if n == 1:
            A[0] += 1
            continue
        facts = factorize(n, spf)
        if all(e == 1 for _, e in facts) and all(p > z for p, _ in facts):
            k = len(facts)
            if k <= 4:
                A[k] += 1

    print(f"\n  A_k counts (factors all > z={z}):")
    for k in range(5):
        print(f"    A_{k}(x,z) = {A[k]}")
    print(f"    Alternating sum A_0 - A_1 + A_2 - A_3 + A_4 = "
          f"{A[0] - A[1] + A[2] - A[3] + A[4]}")
    print(f"    pi(x) - pi(z) = {len(primes) - len(small_primes)}, A_1 = {A[1]}")

    return len(primes), M_x


def explore_transfer_possibilities():
    """
    Explore whether ANY aspect of H-T transfers to pi(x).
    """
    print(f"\n{'='*70}")
    print("EXPLORING TRANSFER POSSIBILITIES")
    print(f"{'='*70}")

    print("""
  POSSIBILITY 1: Express pi(x) as a LINEAR COMBINATION of A_k terms
  ------------------------------------------------------------------
  We know: A_1(x, z) = pi(x) - pi(z)

  Can we express A_1 in terms of other A_k's?
  A_0 = 1 (trivial)
  A_1 = pi(x) - pi(z)  (what we want)
  A_2 = sum_{p > z} (pi(x/p) - pi(p) + 1)  (semiprimes with large factors)
  A_3 = sum_{p > q > z} (pi(x/(pq)) - pi(q) + 1)  (3-almost-primes)
  A_4 = ...

  From inclusion-exclusion or sieve identities:
  A_0 + A_1 + A_2 + A_3 + A_4 = Phi(x, z) = count of n<=x coprime to all p<=z

  So: A_1 = Phi(x, z) - A_0 - A_2 - A_3 - A_4

  This IS the Meissel-Lehmer method! The A_2 term is the "P_2" correction,
  A_3 is "P_3", etc. And computing A_2, A_3, A_4 costs O(x^{2/3}).

  POSSIBILITY 2: Compute pi(x) via M(x) + other summatory functions
  ------------------------------------------------------------------
  pi(x) = sum_{k=1}^{log2(x)} mu(k)/k * Li(x^{1/k}) + corrections

  This relates pi(x) to Li evaluations (easy) and mu(k) values (hard but computable).
  But the EXACT formula requires O(sqrt(x)) terms in the explicit formula.

  If we could compute M(x) in O(x^{3/5}) and somehow convert to pi(x)...
  The conversion itself requires the explicit formula, which costs O(x^{1/2+eps}).
  So the bottleneck would be the CONVERSION, not the M(x) computation.

  POSSIBILITY 3: Signed combination that creates cancellation
  -----------------------------------------------------------
  What if we compute pi(x) - something_easy to get a quantity with cancellation?

  Let S(x) = sum_{n<=x, squarefree} 1 ~ 6x/pi^2  (well-approximated)
  Then pi(x) - S(x) = pi(x) - Q(x) where Q is squarefree count.
  But pi(x) - Q(x) is not a signed sum over multiplicative function.

  What about: pi(x) = sum_{n<=x} Lambda(n)/log(n) + small correction?
  Lambda(n) = log(p) if n = p^k, 0 otherwise.
  This is NOT multiplicative, so Dirichlet convolution techniques don't help.
""")

    print("CONCLUSION:")
    print("  The Helfgott-Thompson technique does NOT transfer to pi(x) because:")
    print("  1. pi(x) counts primes (k=1 almost-primes) -- no signed sum/cancellation")
    print("  2. The H-T trick exploits the alternating sign (-1)^omega in mu")
    print("  3. For pi(x), the Meissel-Lehmer decomposition IS the analogous approach")
    print("     and it costs O(x^{2/3}), not O(x^{3/5})")
    print("  4. Converting M(x) → pi(x) requires the explicit formula, which is")
    print("     at least as expensive as computing pi(x) directly")
    print()
    print("  STATUS: CLOSED. H-T does not transfer to pi(x).")


def main():
    for x in [1000, 10000]:
        analyze_ht_decomposition(x)

    explore_transfer_possibilities()


if __name__ == "__main__":
    main()
