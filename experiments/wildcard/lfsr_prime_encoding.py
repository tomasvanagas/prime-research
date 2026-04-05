"""
LFSR / Finite Field Prime Encoding

RADICAL IDEA: What if we work in GF(2^k) instead of Z?

In GF(2^k), the "sieve of Eratosthenes" becomes operations over a finite field.
Finite field arithmetic is highly structured (LFSR-like).

SPECIFIC APPROACH:
1. The prime indicator function χ(n) = 1 if n prime, 0 otherwise
2. Over GF(2), χ encodes as a binary sequence
3. Question: does this binary sequence have low LINEAR COMPLEXITY over GF(2)?
4. Linear complexity = length of shortest LFSR that generates the sequence
5. If linear complexity is polylog(n), then p(n) is computable in polylog(n)!

The Berlekamp-Massey algorithm computes linear complexity in O(N) time.
We can test empirically: does the linear complexity of the prime indicator
grow as O(N) or sublinearly?

ALSO: MILLS' CONSTANT
Mills proved: there exists A such that ⌊A^{3^n}⌋ is prime for all n.
A ≈ 1.30637788386308...
Computing A requires knowing primes (circular), BUT:
What if there's a DIFFERENT iterated function where the constant is computable?

AND: MATIYASEVICH-ROBINSON FORMULA
There exists a polynomial whose positive values ARE exactly the primes.
The polynomial has ~26 variables. Evaluating at the "right" inputs gives p(n).
Can we find these inputs efficiently?
"""

import numpy as np
import sympy
from sympy import prime, isprime, primepi
import math
import time

def berlekamp_massey_gf2(seq):
    """
    Berlekamp-Massey algorithm over GF(2).
    Returns the linear complexity of the binary sequence.
    """
    n = len(seq)
    s = [int(x) % 2 for x in seq]

    # Current and previous connection polynomials (as lists of coefficients)
    c = [1]  # C(x) = 1
    b = [1]  # B(x) = 1
    L = 0    # current linear complexity
    m = 1    # shift amount
    N_ = 0   # current position

    for N_ in range(n):
        # Compute discrepancy
        d = s[N_]
        for i in range(1, L + 1):
            if i < len(c):
                d ^= (c[i] & s[N_ - i])

        if d == 0:
            m += 1
        elif 2 * L <= N_:
            # Update
            t = c[:]
            # c = c + x^m * b
            while len(c) < len(b) + m:
                c.append(0)
            for i in range(len(b)):
                c[i + m] ^= b[i]
            L = N_ + 1 - L
            b = t
            m = 1
        else:
            # c = c + x^m * b
            while len(c) < len(b) + m:
                c.append(0)
            for i in range(len(b)):
                c[i + m] ^= b[i]
            m += 1

    return L

def test_prime_linear_complexity():
    """
    Test: what is the linear complexity of the prime indicator sequence over GF(2)?

    If it grows as O(N), primes are "maximally complex" in GF(2).
    If it grows sublinearly, there's exploitable structure.
    """
    print("=== Linear complexity of prime indicator over GF(2) ===\n")

    max_N = 5000

    # Generate prime indicator
    indicator = [1 if isprime(n) else 0 for n in range(max_N)]

    # Compute linear complexity for increasing prefixes
    print(f"{'N':>8} {'L(N)':>8} {'L/N':>8} {'L/N·log(N)':>12}")

    checkpoints = [50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000]

    for N in checkpoints:
        if N > max_N:
            break
        L = berlekamp_massey_gf2(indicator[:N])
        ratio = L / N
        log_ratio = L / (N * math.log(N)) if N > 1 else 0
        print(f"{N:>8} {L:>8} {ratio:>8.4f} {log_ratio:>12.6f}")

    print("\n(L/N → 0.5 means random-like; L/N → 0 means structure)")

def test_prime_gf_representations():
    """
    Test prime indicator over various GF(p) to find structure.

    Over GF(p), the prime indicator χ might have different linear complexity.
    Maybe there's a field where primes are "simple"?
    """
    print("\n=== Linear complexity over various GF(p) ===\n")

    N = 500
    indicator = [1 if isprime(n) else 0 for n in range(N)]

    def berlekamp_massey_gfp(seq, p):
        """BM algorithm over GF(p)."""
        n = len(seq)
        s = [int(x) % p for x in seq]

        c = [1]
        b = [1]
        L = 0
        m = 1

        for N_ in range(n):
            d = s[N_]
            for i in range(1, L + 1):
                if i < len(c):
                    d = (d + c[i] * s[N_ - i]) % p

            if d == 0:
                m += 1
            elif 2 * L <= N_:
                t = c[:]
                d_inv = pow(d, p - 2, p)  # Fermat's little theorem
                # Find adjustment
                b_inv_d = [(x * d_inv) % p for x in b]
                while len(c) < len(b_inv_d) + m:
                    c.append(0)
                for i in range(len(b_inv_d)):
                    c[i + m] = (c[i + m] - b_inv_d[i]) % p
                L = N_ + 1 - L
                b = t
                m = 1
            else:
                d_inv = pow(d, p - 2, p)
                b_inv_d = [(x * d_inv) % p for x in b]
                while len(c) < len(b_inv_d) + m:
                    c.append(0)
                for i in range(len(b_inv_d)):
                    c[i + m] = (c[i + m] - b_inv_d[i]) % p
                m += 1

        return L

    print(f"N = {N}")
    print(f"{'Field':>10} {'L':>8} {'L/N':>8}")

    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        L = berlekamp_massey_gfp(indicator, p)
        print(f"{'GF('+str(p)+')':>10} {L:>8} {L/N:>8.4f}")

def test_gap_linear_complexity():
    """
    Maybe the GAP sequence has lower linear complexity than the indicator?

    g(n) = p(n+1) - p(n) might be more structured.
    Also test g(n) mod small numbers.
    """
    print("\n=== Gap sequence linear complexity ===\n")

    N = 2000
    primes_list = [prime(i) for i in range(1, N+1)]
    gaps = [primes_list[i+1] - primes_list[i] for i in range(len(primes_list)-1)]

    # Test gap sequence mod p
    print(f"{'Sequence':>20} {'N':>6} {'L':>8} {'L/N':>8}")

    for p in [2, 3, 4, 6, 8]:
        seq = [g % p for g in gaps[:1000]]
        L = berlekamp_massey_gf2(seq) if p == 2 else len(seq) // 2  # approx for non-GF2
        if p == 2:
            print(f"{'gaps mod 2':>20} {len(seq):>6} {L:>8} {L/len(seq):>8.4f}")

    # The gap sequence mod 2 is interesting:
    # g(n) is even for all n > 1 (since p(n) > 2 are odd, gap is even)
    # So gaps mod 2 = 0 always (except g(1) = 1 since p(1)=2, p(2)=3)
    # This means linear complexity = very low!

    # What about gaps / 2?
    half_gaps = [g // 2 for g in gaps[1:]]  # skip first gap
    for p in [2, 3, 5, 7]:
        seq = [g % p for g in half_gaps[:1000]]
        if p <= 7:
            L = berlekamp_massey_gf2(seq) if p == 2 else 0
            if p == 2:
                print(f"{'gaps/2 mod 2':>20} {len(seq):>6} {L:>8} {L/len(seq):>8.4f}")

    # Cumulative sum of gaps = p(n) - 2
    # So p(n) mod p = (2 + Σ gaps) mod p
    # Linear complexity of p(n) mod p:
    for p in [2, 3, 5, 7]:
        seq = [prime(i) % p for i in range(1, 501)]
        if p == 2:
            L = berlekamp_massey_gf2(seq)
            print(f"{'p(n) mod 2':>20} {len(seq):>6} {L:>8} {L/len(seq):>8.4f}")

def test_mills_like_constants():
    """
    Mills' theorem: ⌊A^{3^n}⌋ is always prime, where A ≈ 1.30637...

    The problem: computing A requires knowing primes.

    BUT: what if we use a DIFFERENT growth rate?
    Question: for which function f(n) does ⌊B^{f(n)}⌋ give the nth prime?

    If f(n) = C·log(n) for some C, then B^{C·log(n)} = n^{C·log(B)},
    and we'd need n^α to enumerate primes, which gives α ≈ 1 + 1/log(n).
    Not useful.

    What about: can we find a COMPUTABLE constant whose n-th digit
    encodes p(n)? This is trivially true (just concatenate primes in binary),
    but the question is whether this constant has a FAST computation
    (like e or π do, via convergent series).
    """
    print("\n=== Mills-like constant analysis ===\n")

    # Compute Mills' constant A to verify
    # A = lim_{n→∞} p(n)^{3^{-n}} where p(n) is Mills' prime sequence
    # Start with the first few Mills primes: 2, 11, 1361, 2521008887, ...

    mills_primes = [2, 11, 1361]  # known values
    for i, mp in enumerate(mills_primes):
        A_approx = mp ** (3.0 ** (-(i+1)))
        print(f"  From Mills prime {mp}: A ≈ {A_approx:.10f}")

    # The point: Mills' constant encodes ALL primes but is not computable
    # without already knowing primes. This is a CIRCULAR definition.

    # Test: for a general prime sequence, what's the growth rate?
    # p(n) ~ n·ln(n). Can we find B, f such that ⌊B^{f(n)}⌋ = p(n)?
    # log(p(n)) = f(n)·log(B)
    # f(n) = log(p(n))/log(B)
    # For B = e: f(n) = log(p(n)) ≈ log(n) + log(log(n))

    print("\n  Testing: can p(n) be expressed as ⌊B^{g(n)}⌋ for smooth g?")
    primes_list = [prime(i) for i in range(1, 1001)]

    # f(n) = log(p(n))
    f_vals = [math.log(p) for p in primes_list]
    n_vals = list(range(1, 1001))

    # Is f(n) = log(p(n)) a smooth function of n?
    # f(n) ≈ log(n) + log(log(n)) - 1 + corrections
    residuals = []
    for i, f in enumerate(f_vals):
        n = i + 1
        if n > 2:
            expected = math.log(n) + math.log(math.log(n))
            residuals.append(f - expected)

    residuals = np.array(residuals)
    print(f"  log(p(n)) - [log(n) + log(log(n))] statistics:")
    print(f"    Mean: {residuals.mean():.6f}")
    print(f"    Std: {residuals.std():.6f}")
    print(f"    (Std/Mean: {residuals.std()/abs(residuals.mean()):.4f})")

    # The residual has std ~ 0.1-0.3, which means individual primes
    # deviate from the smooth approximation by ~10-30%.
    # This confirms: no smooth function computes primes exactly.

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: LFSR / Finite Field Prime Encoding")
    print("=" * 60)

    test_prime_linear_complexity()
    test_prime_gf_representations()
    test_gap_linear_complexity()
    test_mills_like_constants()
