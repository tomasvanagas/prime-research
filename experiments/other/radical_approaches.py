"""
Session 5: Radical new approaches to p(n)

APPROACH 1: "Smooth number gap" method
- Compute R^{-1}(n) to get approximation x₀
- For integers near x₀, test primality using ONLY small primes
- Key insight: most composites near x₀ are divisible by SMALL primes
- After eliminating smooth numbers, very few candidates remain

APPROACH 2: "Prime residue reconstruction"
- p(n) mod small_primes can be determined from π(x; q, a)
- If we could compute π(x; q, a) fast, we could use CRT

APPROACH 3: "Inverse prime counting via inclusion-exclusion on smooth numbers"
- ψ(x, y) = count of x-smooth numbers ≤ x (all prime factors ≤ y)
- π(x) = x - 1 - ψ(x, √x) + ... (inclusion-exclusion)
- Can ψ(x, y) be computed in O(polylog(x))?

APPROACH 4: "Binary representation of p(n)"
- Compute p(n) bit by bit
- For each bit position, determine if it's 0 or 1 using some fast test

APPROACH 5: "Analytic interpolation"
- The function f(n) = p(n) for integer n
- Can we extend to f(z) for complex z and use complex analysis?
- If f(z) is analytic somewhere, Cauchy's integral formula gives exact values

APPROACH 6: "Prime race exploitation"
- Use the precise Chebyshev bias to predict p(n) mod small numbers
- Combine with CRT
"""

from mpmath import mp, mpf, log, exp, li, pi as mpi, sqrt, floor, ceil
from mpmath import nstr, fabs, loggamma, gamma
import sympy
from sympy import prime as sympy_prime, isprime, nextprime, prevprime, primepi
import time
import numpy as np

mp.dps = 50


def approach1_smooth_gap():
    """
    APPROACH 1: Smooth number elimination

    Given approximation x₀ ≈ p(n), the error E means p(n) ∈ [x₀-E, x₀+E].
    In this interval of 2E integers:
    - ~E/2 are even (eliminated)
    - ~E/3 are divisible by 3 (some overlap with even)
    - After sieve by all primes ≤ y, survivors ≈ 2E * ∏(1-1/p) for p≤y

    By Mertens: ∏(1-1/p) ~ e^{-γ}/ln(y)

    So survivors ≈ 2E * e^{-γ} / ln(y)

    For p(10^100): E ≈ 10^34 (using p^0.34 fit)
    survivors after sieve by primes ≤ 10^9: ≈ 2*10^34 / (e^γ * 9*ln(10)) ≈ 10^33

    Still way too many. Even with primes up to 10^20: survivors ≈ 10^32.

    VERDICT: Can't eliminate enough candidates.
    """
    print("=" * 80)
    print("APPROACH 1: SMOOTH NUMBER ELIMINATION")
    print("=" * 80)

    # Test for small n to measure effectiveness
    test_ns = [1000, 10000, 100000]

    for n in test_ns:
        actual = sympy_prime(n)
        # Simulate: error ~ p^0.34
        E = int(actual ** 0.34) + 1
        interval_size = 2 * E

        # Count primes and composites in the interval
        lo = max(2, actual - E)
        hi = actual + E

        # Sieve by small primes
        small_primes = list(sympy.primerange(2, 100))
        survivors = list(range(lo, hi + 1))

        for p in small_primes:
            survivors = [x for x in survivors if x % p != 0 or x == p]

        primes_in_interval = sum(1 for x in range(lo, hi+1) if isprime(x))

        print(f"\nn={n}: p(n)={actual}, E={E}, interval=[{lo},{hi}]")
        print(f"  Interval size: {interval_size}")
        print(f"  After sieve by primes < 100: {len(survivors)} survivors")
        print(f"  Primes in interval: {primes_in_interval}")
        print(f"  Reduction ratio: {len(survivors)/interval_size:.4f}")
        print(f"  Expected by Mertens: {np.exp(-0.5772)/np.log(100):.4f}")


def approach4_bit_by_bit():
    """
    APPROACH 4: Compute p(n) bit by bit.

    p(n) has about log₂(p(n)) ≈ log₂(n * ln(n)) bits.

    For each bit position k (from MSB to LSB):
    - Let x_k = current partial value with bit k set to 0
    - We need: π(x_k + 2^k) - π(x_k) vs expected count

    This requires computing π(x) at each step, which is O(x^{2/3}).
    Total cost: O(log(p(n)) * p(n)^{2/3}) — slightly worse than direct.

    BUT: What if we don't need EXACT π(x) for each bit?
    For the MSB, we need π accurate to ±1, but for intermediate bits,
    we might be able to use approximate π.

    Let's test this idea.
    """
    print("\n" + "=" * 80)
    print("APPROACH 4: BIT-BY-BIT CONSTRUCTION")
    print("=" * 80)

    for target_n in [100, 1000, 10000]:
        actual = sympy_prime(target_n)
        num_bits = actual.bit_length()

        print(f"\nn={target_n}: p(n)={actual} ({num_bits} bits, binary: {bin(actual)})")

        # Build p(n) bit by bit from MSB
        result = 0
        pi_calls = 0

        for bit in range(num_bits - 1, -1, -1):
            # Try setting this bit
            candidate = result | (1 << bit)
            # How many primes ≤ candidate?
            count = primepi(candidate)
            pi_calls += 1

            if count >= target_n:
                # This bit should be 0 or exactly at boundary
                # Check: is primepi(candidate) >= n but primepi(candidate - 2^bit) < n?
                if count == target_n and isprime(candidate):
                    result = candidate
                    # Don't set this bit (keep as is, will be set by later bits)
                    # Actually, we need to think more carefully
                else:
                    pass  # Don't set this bit
            else:
                result = candidate  # Set this bit

        # After all bits, result should be close to actual
        # But we need to find the exact prime

        # Adjust: find the prime at position n
        pi_result = primepi(result)

        print(f"  Bit-by-bit result: {result} (π({result}) = {pi_result})")
        print(f"  Actual: {actual}")
        print(f"  π calls: {pi_calls}")
        print(f"  Match: {result == actual}")

        # Better approach: binary search on x such that π(x) = n
        lo, hi = 2, 2 * actual
        bs_calls = 0
        while lo < hi:
            mid = (lo + hi) // 2
            if primepi(mid) < target_n:
                lo = mid + 1
            else:
                hi = mid
            bs_calls += 1

        print(f"  Binary search result: {lo} (calls: {bs_calls})")
        print(f"  Is prime: {isprime(lo)}")

        # The lo should be p(n) if we search correctly
        # Actually, we want smallest x with π(x) >= n
        # Then we need x to be prime and π(x) = n
        while not isprime(lo):
            lo -= 1
        print(f"  Adjusted: {lo}, correct: {lo == actual}")


def approach5_analytic_continuation():
    """
    APPROACH 5: Analytic continuation of p(n)

    Can we define P(z) for complex z such that P(n) = p(n)?

    One natural extension: P(z) = R^{-1}(z) (inverse Riemann R function)
    This is analytic for Re(z) > some threshold.

    But we know R^{-1}(n) ≠ p(n) exactly.

    Alternative: use the Dirichlet series generating function
    F(s) = Σ p(n) * n^{-s}

    This should converge for Re(s) > 2 (since p(n) ~ n*ln(n)).
    Can we evaluate F at specific points to recover p(n)?

    By Perron's formula: p(N) = (1/2πi) ∫ F(s) * N^s / s ds (contour integral)
    But this requires knowing F(s) already.

    Another idea: Hardy-Littlewood Tauberian theorem.
    """
    print("\n" + "=" * 80)
    print("APPROACH 5: ANALYTIC CONTINUATION")
    print("=" * 80)

    # Compute the Dirichlet series F(s) = Σ p(n)/n^s numerically
    from mpmath import zeta, nsum, inf

    # F(s) for real s > 2
    def prime_dirichlet(s, N=10000):
        """Compute Σ_{n=1}^{N} p(n) * n^{-s}"""
        result = mpf(0)
        for n in range(1, N + 1):
            pn = sympy_prime(n)
            result += mpf(pn) * mpf(n) ** (-mpf(s))
        return result

    # Test: can F(s) be expressed in terms of known functions?
    print("\nDirichlet series F(s) = Σ p(n)/n^s:")
    for s in [3, 4, 5, 10]:
        val = prime_dirichlet(s, 1000)
        print(f"  F({s}) = {nstr(val, 15)}")

    # Check: does F(s) have a simple relationship to ζ(s)?
    # P(s) = Σ p^{-s} (prime zeta function)
    # Our F(s) = Σ p(n) * n^{-s} is DIFFERENT from P(s)

    # Note: p(n) ~ n*ln(n), so F(s) ~ Σ n*ln(n)/n^s = Σ ln(n)/n^{s-1}
    # = -ζ'(s-1) (derivative of Riemann zeta)

    print("\nComparison with -ζ'(s-1):")
    for s in [3, 4, 5, 10]:
        from mpmath import diff
        neg_zeta_prime = -diff(lambda x: zeta(x), s - 1)
        f_val = prime_dirichlet(s, 5000)
        ratio = f_val / neg_zeta_prime
        print(f"  s={s}: F(s)={nstr(f_val,12)}, -ζ'(s-1)={nstr(neg_zeta_prime,12)}, ratio={nstr(ratio,8)}")


def approach6_prime_race():
    """
    APPROACH 6: Use prime race (Chebyshev bias) to predict p(n) mod q.

    For q=4: typically π(x;4,3) > π(x;4,1) (Chebyshev bias)
    For q=3: typically π(x;3,2) > π(x;3,1)

    Quantitative form (Rubinstein-Sarnak):
    δ(x;q,a,b) = (π(x;q,a) - π(x;q,b)) / (√x / ln(x))

    This follows a specific distribution related to zeros of L(s,χ).

    If we could predict δ precisely enough, we'd know π(x;q,a) for each class,
    hence we could determine p(n) mod q.
    """
    print("\n" + "=" * 80)
    print("APPROACH 6: PRIME RACE / CHEBYSHEV BIAS")
    print("=" * 80)

    # Compute prime race statistics
    from collections import defaultdict

    for q in [3, 4, 5, 6, 7, 8, 12]:
        counts = defaultdict(int)
        for p in sympy.primerange(2, 100000):
            counts[p % q] += 1

        total = sum(counts.values())
        print(f"\nPrimes mod {q} (up to 100000):")
        for r in sorted(counts.keys()):
            if counts[r] > 0:
                expected = total / (q - sum(1 for d in range(1, q) if sympy.gcd(d, q) > 1))
                bias = counts[r] / expected if expected > 0 else 0
                from math import gcd
                if gcd(r, q) == 1 or r == 0:
                    print(f"  residue {r}: {counts[r]} ({counts[r]/total*100:.1f}%), bias={bias:.4f}")


def approach7_prime_encoding_entropy():
    """
    APPROACH 7: Information-theoretic analysis of prime encodings.

    The prime gap sequence g(n) = p(n+1) - p(n) has specific entropy.
    If we can find a LOWER entropy encoding, it might be computable.

    Known: g(n) is even for n > 1, median ~ ln(p(n))
    Hardy-Littlewood conjecture gives specific gap probabilities.

    Question: What is the entropy of the gap sequence? How does it compare
    to the information needed to specify p(n)?
    """
    print("\n" + "=" * 80)
    print("APPROACH 7: GAP SEQUENCE ENTROPY ANALYSIS")
    print("=" * 80)

    from collections import Counter
    import math

    # Compute gaps for first N primes
    N = 100000
    primes = list(sympy.primerange(2, sympy_prime(N) + 10))[:N]
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

    # Gap distribution
    gap_counts = Counter(gaps)
    total_gaps = len(gaps)

    print(f"First {N} prime gaps:")
    print(f"  Mean gap: {np.mean(gaps):.3f} (expected: ln({primes[-1]}) = {np.log(primes[-1]):.3f})")
    print(f"  Median gap: {np.median(gaps):.0f}")
    print(f"  Max gap: {max(gaps)}")

    # Shannon entropy
    entropy = 0
    for g, count in gap_counts.items():
        p = count / total_gaps
        if p > 0:
            entropy -= p * math.log2(p)

    print(f"\n  Shannon entropy of gap: {entropy:.4f} bits per gap")
    print(f"  Entropy of uniform on [1, max_gap]: {math.log2(max(gaps)):.4f} bits")
    print(f"  Information to specify p(n): ~log2(p(n)) ≈ {math.log2(primes[-1]):.1f} bits")
    print(f"  Total info in gap sequence: {entropy * total_gaps:.0f} bits")
    print(f"  Direct encoding of p(N): {math.log2(primes[-1]):.1f} bits")

    # Conditional entropy: H(g(n) | g(n-1))
    # Compute transition probabilities
    pair_counts = Counter()
    for i in range(len(gaps)-1):
        pair_counts[(gaps[i], gaps[i+1])] += 1

    cond_entropy = 0
    for (g1, g2), count in pair_counts.items():
        p_pair = count / (total_gaps - 1)
        p_g1 = gap_counts[g1] / total_gaps
        if p_pair > 0 and p_g1 > 0:
            cond_entropy -= p_pair * math.log2(p_pair / p_g1)

    print(f"\n  Conditional entropy H(g_n | g_{{n-1}}): {cond_entropy:.4f} bits")
    print(f"  Reduction from conditioning: {entropy - cond_entropy:.4f} bits ({(entropy-cond_entropy)/entropy*100:.1f}%)")

    # Higher-order conditional: H(g_n | g_{n-1}, g_{n-2})
    triple_counts = Counter()
    for i in range(len(gaps)-2):
        triple_counts[(gaps[i], gaps[i+1], gaps[i+2])] += 1

    cond2_entropy = 0
    for (g1, g2, g3), count in triple_counts.items():
        p_triple = count / (total_gaps - 2)
        p_g1g2 = pair_counts.get((g1, g2), 0) / (total_gaps - 1)
        if p_triple > 0 and p_g1g2 > 0:
            cond2_entropy -= p_triple * math.log2(p_triple / p_g1g2)

    print(f"  H(g_n | g_{{n-1}}, g_{{n-2}}): {cond2_entropy:.4f} bits")
    print(f"  Additional reduction: {cond_entropy - cond2_entropy:.4f} bits")

    # Key question: is the gap sequence COMPRESSIBLE?
    # If H(gap) << log2(max_gap), then there's redundancy we can exploit
    print(f"\n  Compression ratio: {entropy / math.log2(max(gaps)):.4f}")
    print(f"  (1.0 = incompressible, lower = more compressible)")

    # Cumulative gap approach: can we predict sum of gaps?
    cumulative = np.cumsum(gaps)
    # p(n) = 2 + cumulative[n-2] for n >= 2

    # The cumulative sum has much lower relative variance
    # Var(sum) = n * Var(gap) * (1 + 2*sum_autocorrelations)
    gap_var = np.var(gaps)
    print(f"\n  Gap variance: {gap_var:.2f}")
    print(f"  Std of cumulative sum at n={N}: {np.sqrt(N * gap_var):.2f}")
    print(f"  Relative std: {np.sqrt(N * gap_var) / cumulative[-1] * 100:.4f}%")


if __name__ == "__main__":
    approach1_smooth_gap()
    approach4_bit_by_bit()
    approach6_prime_race()
    approach7_prime_encoding_entropy()
    # approach5 is slow, skip for now
