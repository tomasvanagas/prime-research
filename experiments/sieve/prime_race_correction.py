"""
Session 6: Prime Race / Chebyshev Bias Self-Correction

The "prime race" refers to the observation that π(x; 4, 3) > π(x; 4, 1)
for most x (primes ≡ 3 mod 4 tend to lead those ≡ 1 mod 4).

This bias has magnitude O(sqrt(x)/ln(x)) — the SAME order as the
prime counting error!

IDEA: Can we use MULTIPLE prime races (mod 3, mod 4, mod 5, mod 7, ...)
to triangulate the exact error in R^{-1}(n)?

The error in π(x) is: π(x) - li(x) = -Σ_ρ li(x^ρ) + small terms
The error in π(x; q, a) - li(x)/φ(q) involves character sums over
zeros of Dirichlet L-functions.

If we compute MULTIPLE such biases, we get a system of equations
involving DIFFERENT L-function zeros. Maybe the combined system
is easier to solve?
"""

import numpy as np
from sympy import prime, primepi, isprime
from collections import Counter
import time

def compute_prime_race_bias(max_x, q, a_values):
    """
    Compute the bias π(x; q, a) - π(x; q, b) for various residue pairs.
    """
    counts = {a: 0 for a in a_values}
    biases = {}

    for p_val in range(2, max_x + 1):
        if isprime(p_val):
            r = p_val % q
            if r in counts:
                counts[r] += 1

    total = sum(counts.values())
    for a in a_values:
        biases[a] = counts[a] - total / len(a_values)

    return counts, biases

def experiment_multiple_races():
    """
    Use multiple prime races to build a correction system.
    """
    print("="*70)
    print("MULTIPLE PRIME RACE CORRECTIONS")
    print("="*70)

    N = 5000
    primes_list = [int(prime(n)) for n in range(1, N + 1)]

    # For each prime, compute its residue class membership
    moduli = [3, 4, 5, 7, 8, 11, 12]

    # Build residue matrix
    print(f"\nPrime residue biases (first {N} primes):")
    for q in moduli:
        coprime_residues = [r for r in range(1, q) if all(r % p != 0 for p in [2, 3, 5, 7] if p < q and q % p != 0)]
        # Simpler: just count
        counts = Counter(p % q for p in primes_list if p > q)
        total = sum(counts.values())

        print(f"\n  mod {q}:")
        for r in sorted(counts.keys()):
            if counts[r] > 0:
                expected = total / len([rr for rr in range(q) if np.gcd(rr, q) == 1]) if any(np.gcd(rr, q) == 1 for rr in range(q)) else total
                bias = counts[r] - total / max(1, sum(1 for rr in range(q) if np.gcd(rr, q) == 1))
                print(f"    r={r}: count={counts[r]}, bias={bias:+.1f}")

    # Key experiment: Can residue class biases predict the error δ(n)?
    # For each n, the residue p(n) mod q is either determined by the formula
    # or provides information about which candidate prime is correct.

    print(f"\n\nCan residue biases correct the approximation?")
    from numpy import gcd

    correct_base = 0
    correct_residue = 0

    for n in range(10, N + 10):
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n)

        # Base: round to nearest prime
        x = int(round(approx))
        from sympy import nextprime, prevprime
        cands = set()
        try:
            cands.add(int(nextprime(x - 1)))
            cands.add(int(prevprime(x + 1)))
            for offset in range(-30, 31, 2):
                c = x + offset
                if c > 1 and isprime(c):
                    cands.add(c)
        except:
            pass

        if p_n in cands:
            # Among candidates, can residue information help?
            # Score each candidate by how well it matches expected biases
            best_cand = min(cands, key=lambda c: abs(c - approx))
            if best_cand == p_n:
                correct_base += 1

            # Try scoring by Chebyshev bias:
            # primes ≡ 3 (mod 4) should slightly outnumber ≡ 1 (mod 4)
            # So if candidate ≡ 3 (mod 4), it's slightly more likely
            scores = {}
            for c in cands:
                score = -abs(c - approx)  # Distance penalty
                # Chebyshev bias bonuses (very small)
                if c % 4 == 3:
                    score += 0.1
                if c % 3 == 2:
                    score += 0.05
                scores[c] = score

            best_biased = max(scores, key=scores.get)
            if best_biased == p_n:
                correct_residue += 1

    total = N
    print(f"  Base (nearest to approx): {correct_base}/{total} ({100*correct_base/total:.1f}%)")
    print(f"  With residue bias:        {correct_residue}/{total} ({100*correct_residue/total:.1f}%)")
    print(f"  Improvement: {correct_residue - correct_base} ({100*(correct_residue-correct_base)/max(1,total):.1f}%)")

    # The Chebyshev bias is O(1/ln(p)) per prime, which is tiny
    # It can't meaningfully help distinguish between nearby candidates
    print(f"\n  CONCLUSION: Chebyshev bias is O(1/ln(p)) ≈ {1/np.log(primes_list[N//2]):.4f}")
    print(f"  This is far too weak to distinguish candidates separated by O(1)")

def main():
    print("Session 6: Prime Race / Chebyshev Bias")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()
    experiment_multiple_races()
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
