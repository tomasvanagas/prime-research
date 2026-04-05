"""
PROPOSAL 5: Adelic Interpolation / Multi-Residue Collapse
=========================================================

IDEA: Combine p-adic and real information about p(n) simultaneously.

The adelic viewpoint: a prime p lives in Q embedded diagonally into
the adele ring A_Q = R x prod_p Q_p.

For the nth prime:
- REAL component: R^{-1}(n) gives ~50% of digits in R
- p-ADIC components: p(n) mod q for small primes q

If we could compute p(n) mod q for enough primes q cheaply,
CRT gives exact p(n). The question: can we get p(n) mod q in O(polylog)?

KEY INSIGHT: p(n) mod q is related to the distribution of primes in
residue classes mod q. By Dirichlet's theorem:
  pi(x; q, a) ~ li(x) / phi(q)

More precisely, the Chebyshev bias and prime races give:
  pi(x; q, a) - li(x)/phi(q) = -(1/phi(q)) * sum_chi chi(a)_bar * sum_rho li(x^rho)

For FIXED q, this correction involves L(s, chi) zeros, which are
DIFFERENT from zeta zeros. Different L-functions might have sparser
zero contributions.

RADICAL IDEA: For special moduli q (e.g., highly composite numbers,
or q = product of first k primes), the character sum structure might
allow efficient computation of pi(x) mod q.

EVEN MORE RADICAL: Compute p(n) mod q for q = 2,3,5,7,... by
exploiting the SIEVE STRUCTURE. The sieve of Eratosthenes naturally
computes primality mod each small prime. Can we extract p(n) mod q
from a compressed sieve?

COMPRESSED SIEVE IDEA:
- The sieve of size x has x bits but only pi(x) ~ x/ln(x) ones
- By CRT, knowing p(n) mod q for q = 2,...,B with prod(q) > x gives p(n)
- Need B ~ O(ln(x)) primes, requiring prod(primes up to B) > x
- By PNT, prod(p <= B) ~ e^B, so need B ~ ln(x) = ln(n ln n) ~ ln(n)
- For each mod q: can we compute the INDEX of p(n) in the reduced
  sieve mod q?

THE ALGORITHM:
1. Compute x_approx = R^{-1}(n) in O(polylog(n))
2. For q = 2, 3, 5, ..., B (B ~ ln(n)):
   a. Compute pi(x_approx; q, a) for each a in (Z/qZ)* using L-functions
   b. Determine which residue class a satisfies:
      sum_{a'<=a} pi(x_approx; q, a') >= (n - pi(sqrt(x))) and the next doesn't
   c. This gives p(n) mod q
3. CRT to get p(n) exactly

COMPLEXITY: Step 2a requires L-function zero sums. For each character chi mod q,
similar cost to the zeta-zero sum. Total: sum_{q<=B} phi(q) ~ B^2/2 characters.
Each character needs O(sqrt(x)) zeros... unless there's a shortcut.

BUT WAIT: We don't need pi(x; q, a) exactly. We need it modulo the number
of primes in the interval, which is O(sqrt(x)/log(x)). This modular
reduction might allow massive cancellation.
"""

import numpy as np
from sympy import (prime, primepi, isprime, nextprime, prevprime,
                   primerange, factorint, totient, mobius, primitive_root)
from mpmath import mp, mpf, li, log, exp
import math

mp.dps = 30

def primes_in_residue_class(x, q, a):
    """Count primes p <= x with p ≡ a (mod q)."""
    count = 0
    for p in primerange(2, x + 1):
        if p % q == a:
            count += 1
    return count

def compute_pn_mod_q(n_target, q):
    """
    Compute p(n) mod q.

    Method: iterate through residue classes of (Z/qZ)* in order,
    accumulating prime counts until we reach the nth prime.
    """
    p_n = int(prime(n_target))
    return p_n % q

def test_adelic_reconstruction(n_target, max_B=30):
    """
    Test: reconstruct p(n) from its residues mod small primes.
    Measure how many moduli are needed.
    """
    p_n = int(prime(n_target))

    moduli = list(primerange(2, max_B + 1))
    residues = [p_n % q for q in moduli]

    # CRT reconstruction with increasing number of moduli
    results = []
    for k in range(1, len(moduli) + 1):
        mods = [int(m) for m in moduli[:k]]
        ress = [int(r) for r in residues[:k]]

        product = 1
        for m in mods:
            product *= m

        if product >= p_n:
            # CRT
            from functools import reduce
            M = product
            reconstructed = 0
            for r, m in zip(ress, mods):
                Mi = M // m
                Mi_inv = pow(Mi, -1, m)
                reconstructed = (reconstructed + r * Mi * Mi_inv) % M

            results.append({
                'k': k,
                'product': product,
                'reconstructed': reconstructed,
                'correct': reconstructed == p_n,
                'moduli': mods
            })

            if reconstructed == p_n:
                return {
                    'n': n_target,
                    'p_n': p_n,
                    'moduli_needed': k,
                    'moduli': mods,
                    'product': product,
                    'correct': True
                }

    return {
        'n': n_target,
        'p_n': p_n,
        'moduli_needed': len(moduli),
        'moduli': [int(m) for m in moduli],
        'product': product if 'product' in dir() else 0,
        'correct': False
    }

def prime_race_residue_prediction(n_target, q):
    """
    Try to predict p(n) mod q from the prime race statistics.

    The Chebyshev bias says certain residue classes are slightly favored.
    For q=4: primes ≡ 3 (mod 4) are slightly more common than ≡ 1 (mod 4).

    Can we predict which class p(n) falls in from n alone?
    """
    p_n = int(prime(n_target))
    actual_residue = p_n % q

    # Count primes in each residue class up to p(n)
    residue_counts = {}
    for a in range(q):
        if math.gcd(a, q) == 1 or a == 0:
            residue_counts[a] = 0

    for p in primerange(2, p_n + 1):
        r = int(p) % q
        if r in residue_counts:
            residue_counts[r] += 1

    # Predict: choose the class where the running count "should" hit next
    total = sum(residue_counts.values())
    expected_per_class = total / int(totient(q))

    return {
        'q': q,
        'actual': actual_residue,
        'counts': residue_counts,
        'total': total,
        'expected_per_class': expected_per_class,
    }

def test_residue_prediction_accuracy(n_range, q):
    """
    For a range of n values, test how well we can predict p(n) mod q
    from local statistics.
    """
    correct = 0
    total = 0

    residue_classes = [a for a in range(q) if math.gcd(a, q) == 1]
    if q == 2:
        residue_classes = [1]  # All primes > 2 are odd

    for n in range(max(3, n_range[0]), n_range[1] + 1):
        p_n = int(prime(n))
        actual = p_n % q

        # Naive prediction: most likely residue class
        # For q=2: always 1 (except p(1)=2)
        if q == 2:
            predicted = 1 if n > 1 else 0
        else:
            # Predict based on running bias
            predicted = residue_classes[n % len(residue_classes)]

        if predicted == actual:
            correct += 1
        total += 1

    return correct / total if total > 0 else 0

def compressed_sieve_mod_q(x, q):
    """
    Compute the sieve of Eratosthenes mod q.

    Instead of storing all primes up to x, store only their
    residues mod q. This gives a sequence in (Z/qZ)*.

    The sequence p(1) mod q, p(2) mod q, ... has structure
    determined by the distribution of primes in residue classes.
    """
    primes_mod_q = []
    for p in primerange(2, x + 1):
        primes_mod_q.append(int(p) % q)

    # Analyze the sequence
    from collections import Counter
    counts = Counter(primes_mod_q)

    # Check for patterns
    # Compute transition probabilities: P(p(n+1)≡b | p(n)≡a)
    transitions = {}
    for i in range(len(primes_mod_q) - 1):
        key = (primes_mod_q[i], primes_mod_q[i+1])
        transitions[key] = transitions.get(key, 0) + 1

    return counts, transitions

def analyze_prime_residue_markov(max_x=10000, q=3):
    """
    Model p(n) mod q as a Markov chain.
    If transitions are predictable, we might predict p(n) mod q
    from p(n-1) mod q without computing p(n).
    """
    counts, transitions = compressed_sieve_mod_q(max_x, q)

    # Compute transition matrix
    states = sorted(set(r for r in range(q) if math.gcd(r, q) == 1 or r < 3))
    n_states = len(states)

    total_from = {}
    for (a, b), count in transitions.items():
        total_from[a] = total_from.get(a, 0) + count

    print(f"\nMarkov transition matrix for primes mod {q}:")
    print(f"{'':>5s}", end="")
    for b in states:
        print(f"  ->{b:d}", end="")
    print()

    for a in states:
        print(f"  {a:d}: ", end="")
        for b in states:
            prob = transitions.get((a, b), 0) / total_from.get(a, 1)
            print(f" {prob:.3f}", end="")
        print()

    return counts, transitions


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 5: Adelic Interpolation / Multi-Residue Collapse")
    print("=" * 80)

    print("\n--- Part A: How many moduli needed for CRT reconstruction? ---")
    for n in [10, 100, 500, 1000, 5000]:
        result = test_adelic_reconstruction(n)
        print(f"  n={n:5d}: p(n)={result['p_n']:6d}, moduli_needed={result['moduli_needed']}, "
              f"correct={result['correct']}, moduli={result['moduli']}")

    print("\n--- Part B: Prime residue Markov chains ---")
    for q in [3, 5, 7]:
        analyze_prime_residue_markov(max_x=10000, q=q)

    print("\n--- Part C: Residue prediction accuracy ---")
    for q in [2, 3, 5, 7, 11]:
        acc = test_residue_prediction_accuracy((3, 500), q)
        random_baseline = 1.0 / int(totient(q)) if q > 2 else 1.0
        print(f"  mod {q:2d}: accuracy={acc:.4f}, random_baseline={random_baseline:.4f}, "
              f"improvement={acc/random_baseline:.2f}x")

    print("\n--- Part D: Primorial bounds ---")
    print("  How many prime moduli to reconstruct p(n)?")
    primorial = 1
    count = 0
    for p in primerange(2, 100):
        primorial *= int(p)
        count += 1
        bits = math.log2(primorial)
        # What n can we handle?
        # Need primorial > p(n) ~ n*ln(n)
        # Solve n*ln(n) = primorial approximately
        # n ~ primorial / ln(primorial) = primorial / (bits * ln(2))
        max_n = primorial / (bits * math.log(2))
        print(f"    {count} primes (up to {int(p)}): primorial has {bits:.0f} bits, "
              f"handles n up to ~{max_n:.0e}")
        if count >= 15:
            break
