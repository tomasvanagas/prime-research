#!/usr/bin/env python3
"""
Experiment: Algebraic/Arithmetic shortcuts for p(n)

RADICAL IDEAS that bypass analytic number theory entirely:

1. WILSON'S THEOREM APPROACH:
   (p-1)! ≡ -1 (mod p) iff p is prime.
   Can we use modular arithmetic to "detect" primes without trial division?

   χ_P(n) = ((n-1)! + 1) mod n  [0 if prime, nonzero otherwise... almost]

   If we could compute n! mod m fast, we could get prime indicators fast.
   n! mod p is computable in O(√p polylog p) via algorithms of
   Costa, Gerbicz, Harvey (2014).

2. BERNOULLI NUMBER / ZETA VALUE APPROACH:
   ζ(-n) = -B_{n+1}/(n+1), relating zeta to Bernoulli numbers.
   Bernoulli numbers satisfy recurrences. Can we extract prime info from them?

3. CYCLOTOMIC POLYNOMIAL APPROACH:
   The nth cyclotomic polynomial Φ_n(x) has degree φ(n).
   Φ_p(x) = 1 + x + x² + ... + x^{p-1} for prime p.
   Can we use cyclotomic structure to count/locate primes?

4. POLYNOMIAL EVALUATION MOD n:
   For prime p, Z/pZ is a field. Various polynomial identities hold
   only when p is prime. Can we batch-evaluate these?

5. MATRIX PERMANENT / DETERMINANT:
   The number of derangements of multiples relates to the sieve.
   Can we express π(x) as a determinant (computable in O(n^3))?

6. THE FACTORIAL PRIME COUNTING FORMULA:
   π(n) = -1 + Σ_{j=2}^{n} floor((j-1)! - j*floor((j-1)!/j))/(j-1)
   This is exact but O(n) with expensive factorials.
   Can we speed it up?
"""

import numpy as np
from sympy import (primerange, prime, primepi, isprime, factorint,
                   bernoulli, factorial, cyclotomic_poly, Symbol)
from sympy.ntheory import totient
import time

def wilson_based_counting():
    """
    Test: Use Wilson's theorem variants for prime counting.

    Key: χ(n) = 1 - ((n-1)! + 1) mod n  ... but this is 1 only for primes

    Can we compute Σ_{k=2}^{x} χ(k) = π(x) without computing each χ(k)?

    χ(k) involves (k-1)! mod k. By Wilson, this is -1 for prime k, and
    for composite k, (k-1)! ≡ 0 (mod k) [for k ≥ 5, since k = a*b with a,b < k].
    Exception: k=4: 3! = 6 ≡ 2 (mod 4).
    """
    print("=" * 70)
    print("WILSON'S THEOREM APPROACH")
    print("=" * 70)

    # Verify Wilson's theorem
    print("\nVerification:")
    for n in range(2, 30):
        fm = 1
        for k in range(1, n):
            fm = (fm * k) % n
        is_p = (fm == n - 1)
        print(f"  n={n:2d}: ({n-1})! mod {n} = {fm}, "
              f"Wilson says prime: {is_p}, actually prime: {isprime(n)}")

    # Key question: can we compute Σ (k-1)! mod k for k=2..x efficiently?
    #
    # Note: (k)! mod (k+1) = k * (k-1)! mod (k+1)
    # So if we know (k-1)! mod k, can we get (k)! mod (k+1)?
    #
    # (k)! mod (k+1) = k * (k-1)! mod (k+1)
    # But (k-1)! mod k ≠ (k-1)! mod (k+1) in general.

    # The factorial grows too fast for modular tricks to help across different moduli.

    # Alternative: compute (k-1)! mod k using the product formula
    # (k-1)! = 1 * 2 * ... * (k-1)
    # This is O(k) per value of k, O(x²) total. Can we batch?

    # IDEA: Use the Chinese Remainder Theorem.
    # (k-1)! mod k depends on the prime factorization of k.
    # For k = p (prime): (p-1)! ≡ -1 (mod p) [Wilson]
    # For k = p^a: (k-1)! ≡ 0 (mod k) if a ≥ 2 and p^a ≥ 5
    # For k = p*q (distinct primes): (k-1)! ≡ 0 (mod k)

    # So: (k-1)! mod k = { k-1 if k prime, 2 if k=4, 0 if composite k≥5 }
    # This is EXACTLY the primality indicator!
    # But computing it requires factoring k or direct computation.

    print("\n\nConclusion: Wilson's theorem is a RESTATEMENT of primality,")
    print("not a shortcut. Computing (k-1)! mod k is as hard as primality testing.")


def determinant_sieve():
    """
    IDEA: Express the sieve as a determinant.

    Define matrix A where A[i][j] = 1 if prime_j divides i, for i=1..x, j=1..pi(y).
    The number of rows with all zeros = x - rank(A) ... NO, this is wrong.

    Better: Use inclusion-exclusion as a permanent/determinant.

    S(x, {p_1,...,p_k}) = Σ_{S ⊆ {1,...,k}} (-1)^|S| floor(x / prod_{j∈S} p_j)

    This looks like a permanent computation, but with floor functions.

    Alternative: Can we write S(x, {p_1,...,p_k}) as det(I - M) for some matrix M?

    The Euler product: prod(1 - 1/p) is a product of (1 - 1/p) factors.
    Its "matrix analog": det(I - B) where B is block-diagonal with blocks 1/p_j.

    But we need the EXACT count, not just the product.
    """
    print("\n" + "=" * 70)
    print("DETERMINANT SIEVE APPROACH")
    print("=" * 70)

    from itertools import combinations

    for x, y in [(100, 10), (500, 23), (1000, 32)]:
        primes = list(primerange(2, y+1))
        k = len(primes)

        # True sieve count
        true_count = sum(1 for n in range(1, x+1)
                        if all(n % p != 0 for p in primes))

        # Can we express this as a determinant?
        #
        # Consider the matrix M of size k × k where:
        # M[i][j] = floor(x / (p_i * p_j)) for i ≠ j, M[i][i] = floor(x / p_i)
        #
        # Then det(I - M/x) ≈ prod(1 - 1/p_i) ... the Euler product.
        # But this is just the DENSITY, not the exact count.

        # BETTER IDEA: Smith normal form.
        # The GCD matrix G[i][j] = gcd(i,j) has known determinant = prod φ(k).
        # Can we use GCD matrices to count coprimes?

        # The number of n ≤ x coprime to m = Σ_{d|m} μ(d) floor(x/d)
        # For m = primorial, this is a signed sum of floor functions.

        # MATRIX REPRESENTATION:
        # Define vector v = (floor(x/d))_{d | m} and vector μ = (μ(d))_{d | m}
        # Then S(x, primes) = v · μ = dot product.
        # This is O(2^k) for the dot product, not a determinant.

        # Can we compute v · μ in O(poly(k)) time?
        # v · μ = Σ_{d|m} μ(d) floor(x/d)
        # = x * Σ_{d|m} μ(d)/d - Σ_{d|m} μ(d) {x/d}
        # = x * φ(m)/m - Σ_{d|m} μ(d) {x/d}

        # The first term is O(k) to compute. The second is the hard part.
        # There are 2^k divisors of m, so the second sum has 2^k terms.

        # KEY QUESTION: Is Σ_{d|m} μ(d) {x/d} computable in poly(k) time?

        # For m = p1*p2*...*pk (squarefree):
        # Divisors d = product of a subset of primes.
        # {x/d} = x/d - floor(x/d)

        # Σ μ(d){x/d} = Σ_{S⊆[k]} (-1)^|S| {x/∏_{i∈S} p_i}

        # This is a SIGNED SUM of FRACTIONAL PARTS. The fractional parts
        # are pseudo-random (equidistributed by Weyl's theorem for irrational x/d).

        # For a sum of 2^k pseudo-random signs, the expected magnitude is 2^{k/2}.
        # This is << 2^k terms but >> 1.

        # Can we compute this via Fourier analysis on the group (Z/2)^k?

        frac_sum = 0
        term_count = 0
        for size in range(k+1):
            for S in combinations(range(k), size):
                d = 1
                for i in S:
                    d *= primes[i]
                mu_d = (-1)**size
                frac = x/d - x//d
                frac_sum += mu_d * frac
                term_count += 1

        euler_approx = x * np.prod([1 - 1/p for p in primes])

        print(f"\nx={x}, y={y}, k={k}")
        print(f"  True count: {true_count}")
        print(f"  Euler approx: {euler_approx:.2f}")
        print(f"  Fractional sum: {frac_sum:.6f}")
        print(f"  |frac_sum| / 2^{k//2}: {abs(frac_sum) / 2**(k/2):.4f}")
        print(f"  Terms computed: {term_count} = 2^{k}")

    print("\n  CONCLUSION: The fractional sum has magnitude O(2^{k/2}),")
    print("  consistent with pseudo-random cancellation.")
    print("  No obvious way to compute it in poly(k) time.")


def cyclotomic_approach():
    """
    IDEA: Cyclotomic polynomials encode prime structure.

    Φ_n(1) = { p if n = p^k for prime p
              { 1 otherwise

    So: Π_{n=1}^{x} Φ_n(1) = Π_{prime powers p^k ≤ x} p

    And: log(Π Φ_n(1)) = Σ_{p^k ≤ x} log(p) = ψ(x) (Chebyshev's function)

    But this requires evaluating x cyclotomic polynomials at 1, which is O(x).

    What if we work at a different point? Φ_n(z) for z ≠ 1?
    """
    print("\n" + "=" * 70)
    print("CYCLOTOMIC POLYNOMIAL APPROACH")
    print("=" * 70)

    # Verify: Φ_n(1)
    x = Symbol('x')
    print("\nΦ_n(1) values:")
    for n in range(1, 31):
        phi_n_1 = cyclotomic_poly(n, 1)
        print(f"  n={n:2d}: Φ_n(1) = {phi_n_1}", end="")
        if phi_n_1 > 1:
            print(f"  ← {phi_n_1} is prime, n is a power of {phi_n_1}")
        else:
            print()

    # Product formula: Π_{d|n} Φ_d(x) = x^n - 1
    # So: Π_{d|n} Φ_d(1) = 0 for all n > 1.
    # Actually Π_{d|n} Φ_d(x) = x^n - 1, so at x=1: Π_{d|n} Φ_d(1) = 0.
    # This means Φ_1(1) = 0 (since Φ_1(x) = x-1).
    # For n > 1: n = Π_{d|n, d>1} Φ_d(1) ... wait, that's not right.

    # Actually: x^n - 1 = Π_{d|n} Φ_d(x)
    # Taking derivative and evaluating at x=1: n = Σ_{d|n} Π_{d'|n, d'≠d} Φ_{d'}(1)
    # This involves Φ_1(1) = 0, making things complicated.

    # Let's try x = 2 instead.
    print("\nΦ_n(2) values and their factorizations:")
    product = 1
    for n in range(1, 21):
        phi_n_2 = cyclotomic_poly(n, 2)
        product *= phi_n_2
        # Π_{d|n} Φ_d(2) = 2^n - 1
        check = 2**n - 1
        print(f"  n={n:2d}: Φ_n(2) = {str(phi_n_2):>8s}, "
              f"Π_{{d|n}} Φ_d(2) = {check}")

    # Mersenne numbers: 2^n - 1 = Π_{d|n} Φ_d(2)
    # If n is prime: 2^n - 1 = Φ_1(2) * Φ_n(2) = 1 * Φ_n(2) = Φ_n(2)
    # So Φ_p(2) = 2^p - 1 (Mersenne number)

    # This doesn't directly help with prime counting.

    # IDEA: Is there a generating function using cyclotomic polynomials
    # that encodes π(x)?

    # Π_{p ≤ x} (1 - z/p) = ... involves primes directly.
    # But computing this product for specific z requires knowing all primes ≤ x.

    print("\n  CONCLUSION: Cyclotomic polynomials encode prime structure")
    print("  but don't provide a shortcut for counting primes.")


def fast_factorial_mod():
    """
    Test: How fast can we compute n! mod p using modern algorithms?

    The best known: O(p^{1/2} polylog p) by Harvey (2014).
    This was later improved to O(p^{1/2} (log p)^{1+o(1)}).

    If we could compute k! mod p for all k ≤ x simultaneously,
    we could identify primes via Wilson's theorem.

    But this is a per-prime computation, not a batch operation.

    NOVEL ANGLE: Can we use the structure of n! mod (n+1) across many n?

    n! mod (n+1) has a specific relationship to the factorization of n+1.
    For prime n+1: n! ≡ -1 (mod n+1) [Wilson]
    For n+1 = p*q: n! ≡ 0 (mod n+1) [since both p and q appear in n!]
    """
    print("\n" + "=" * 70)
    print("FAST FACTORIAL MOD APPROACH")
    print("=" * 70)

    # Compute n! mod (n+1) incrementally
    # This is O(n) total since we maintain a running factorial

    print("\nIncremental factorial mod computation:")
    fact_mod = 1  # 1! mod 2
    primes_found = []

    x = 1000
    t0 = time.time()
    for n in range(1, x):
        # We have (n-1)! mod n from previous iteration... but moduli change!
        # Can't reuse directly.
        pass

    # Actually, the issue is that n! mod (n+1) and (n+1)! mod (n+2)
    # use DIFFERENT moduli. So we can't maintain a running product.

    # We need: (n-1)! mod n for each n.
    # (n-1)! = (n-2)! * (n-1)
    # But (n-2)! mod (n-1) ≠ (n-2)! mod n.

    # To get (n-1)! mod n from (n-2)! mod (n-1):
    # We need (n-2)! mod n, which we don't have.

    # RECURSIVE APPROACH: Maintain (k)! mod n for all n simultaneously?
    # That requires O(x) storage and O(x) work per k, giving O(x²) total.
    # No better than brute force.

    # BATCH APPROACH: Compute n! mod n for all n ≤ x using FFT?
    # n! mod n = 0 for all n ≥ 2. Not useful.
    # (n-1)! mod n: this requires modular arithmetic with changing modulus.

    # KEY INSIGHT: For the Wilson sieve, we need (n-1)! mod n for each n.
    # This is essentially asking "is n prime?" for each n.
    # Any batch algorithm for this IS a prime sieve.

    # The Sieve of Eratosthenes is O(x log log x), which is faster than
    # any factorial-based approach could be.

    print("  Computing (n-1)! mod n for n = 2 to 100 via direct computation:")
    for n in range(2, 101):
        fm = 1
        for k in range(1, n):
            fm = (fm * k) % n
        if fm == n - 1:
            primes_found.append(n)

    print(f"  Primes found: {primes_found[:20]}...")
    print(f"  Total: {len(primes_found)}")
    print(f"  Expected: {primepi(100)}")

    print("\n  CONCLUSION: Wilson's theorem gives an O(x²) prime sieve,")
    print("  much worse than Eratosthenes' O(x log log x).")
    print("  No shortcut via factorial arithmetic.")


def matrix_power_approach():
    """
    IDEA: Express the prime indicator as a "state" that evolves
    via matrix multiplication, allowing fast exponentiation.

    For Fibonacci: F_n = M^n * F_0 where M is 2×2.

    For primes: Is there a matrix M of polynomial size such that
    some entry of M^n relates to the nth prime?

    If M were d×d, then M^n costs O(d^3 log n) via repeated squaring.
    For this to be polylog in n, we need d = O(polylog n) AND
    the entries to stay manageable.

    But primes don't satisfy a linear recurrence (the gap sequence
    has full entropy), so no fixed-size matrix M works.

    UNLESS: we work over a ring where the "state" has more structure.
    E.g., matrices over Z/mZ for carefully chosen m.
    """
    print("\n" + "=" * 70)
    print("MATRIX POWER APPROACH")
    print("=" * 70)

    # Test: Does the prime-counting function π(n) satisfy
    # an approximate linear recurrence?

    # π(n) - π(n-1) = 1 if n is prime, 0 otherwise.
    # So: π(n) = π(n-1) + χ_P(n)
    # This is a recurrence, but χ_P(n) is the unknown.

    # What if we work with the VECTOR (π(n), n, 1)?
    # π(n) = π(n-1) + χ_P(n)
    # n = (n-1) + 1
    # 1 = 1

    # Matrix form: [π(n), n, 1]^T = M_n * [π(n-1), n-1, 1]^T
    # where M_n = [[1, 0, χ_P(n)], [0, 1, 1], [0, 0, 1]]

    # The product M_x * M_{x-1} * ... * M_2 applied to [0, 1, 1]^T
    # gives [π(x), x, 1]^T.

    # But M_n depends on whether n is prime, which is what we want to know!

    # WORKAROUND: What if we DON'T condition on primality,
    # but instead carry ALL possible states?

    # State: (count, position) in {0,1,...,k} × {2,...,x}
    # Transition: at position n, if prime, increment count; else don't.
    # This is a 2-state system at each step (prime or not).

    # After x steps, we've created a binary tree of 2^x paths.
    # The "correct" path has π(x) primes. But we can't select it
    # without knowing which numbers are prime.

    # This is inherently a non-autonomous dynamical system
    # (the transition depends on the current position, not just state).

    # For AUTONOMOUS systems (M constant), fast exponentiation works.
    # For non-autonomous (M varies), the product has no shortcut in general.

    print("Testing if the prime indicator sequence has approximate periodicity...")

    # Check: for various moduli m, does the sequence χ_P(n) mod m
    # (which is just χ_P(n) since it's 0 or 1) have patterns?

    # More interesting: does the sequence of gaps g(n) = p(n+1) - p(n)
    # have periodicity mod m?

    primes_list = list(primerange(2, 10000))
    gaps = [primes_list[i+1] - primes_list[i] for i in range(len(primes_list)-1)]

    # Autocorrelation of gaps
    gaps_arr = np.array(gaps, dtype=float)
    gaps_centered = gaps_arr - np.mean(gaps_arr)

    n = len(gaps_centered)
    autocorr = np.correlate(gaps_centered[:min(n, 2000)],
                            gaps_centered[:min(n, 2000)], mode='full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]

    print(f"\nPrime gap autocorrelation (first {len(primes_list)} primes):")
    for lag in [1, 2, 3, 5, 10, 30, 100]:
        if lag < len(autocorr):
            print(f"  lag {lag:4d}: {autocorr[lag]:.4f}")

    # Check for periodicity in gaps mod small numbers
    print(f"\nPrime gap distribution mod small numbers:")
    for m in [2, 3, 6, 12, 30]:
        counts = np.zeros(m, dtype=int)
        for g in gaps:
            counts[g % m] += 1
        print(f"  mod {m:2d}: {dict(enumerate(counts))}")

    # KEY FINDING: Gaps mod 6 are always 0 (for twin primes: g=2,4)
    # or multiples of 6 ± patterns. This is because primes > 3 are ≡ 1 or 5 mod 6.

    print("\n  CONCLUSION: Prime gaps have local correlations (Hardy-Littlewood)")
    print("  but no periodic structure suitable for matrix exponentiation.")
    print("  The gap sequence has positive entropy, ruling out fixed-size matrices.")


def quadratic_form_approach():
    """
    IDEA: Represent prime counting as optimizing a quadratic form.

    The Selberg sieve: minimize Σ_{n≤x} (Σ_{d|n} λ_d)² subject to λ_1 = 1.
    The optimal λ values give upper bounds for π(x).

    What if we set up an EXACT formulation as a quadratic program?

    Let v ∈ {0,1}^x where v_n = 1 iff n is prime.
    Then π(x) = 1^T v.

    Constraints: for each composite n = ab (a,b > 1), v_n = 0.
    These are LINEAR constraints.

    So: π(x) = max 1^T v subject to v ∈ {0,1}^x, v_n = 0 for composite n.

    This is an integer linear program, but the "composite" constraints
    require knowing which n are composite -- circular again.

    UNLESS we encode compositeness algebraically:
    v_n * v_m * (v_{nm} - 1) ≤ 0 for all n,m with nm ≤ x?
    No, this says "if both n and m are in the set, nm must be too" --
    that's the opposite of what we want.

    We want: for n > 1, if n = ab with a,b > 1, then v_n = 0.
    Equivalently: v_n = 1 implies n has no factorization n = ab with 1 < a ≤ b < n.

    In constraint form: v_n ≤ 1 - v_a * v_b for all a*b = n, 1 < a ≤ b.
    Wait, this is wrong too. v_a and v_b might not be 1.

    Actually: v_n * (Σ_{a*b=n, 1<a≤b} 1) should be 0 if n is composite.
    The number of ways to write n = a*b is related to d(n)-2 (number of divisors minus 2).

    Simpler: v_n * Σ_{d|n, 1<d<n} 1 = 0 for prime v_n.
    This says: if v_n = 1, then n has no divisors other than 1 and n.
    Which IS primality.

    This formulation has x variables and Σ_{n=2}^{x} d(n) ≈ x log x constraints.
    Solving integer programming is NP-hard in general, but this specific
    structure might be exploitable.
    """
    print("\n" + "=" * 70)
    print("QUADRATIC FORM / OPTIMIZATION APPROACH")
    print("=" * 70)

    # Test: Set up the LP relaxation for small x
    # Relax v ∈ {0,1} to v ∈ [0,1]

    x_test = 50

    # Without using LP solver, test the Selberg sieve as a quadratic form
    from itertools import combinations

    primes_up_to_x = list(primerange(2, x_test + 1))
    actual_pi = len(primes_up_to_x)

    print(f"\nFor x = {x_test}:")
    print(f"  Actual π(x) = {actual_pi}")

    # Selberg upper bound: π(x) ≤ 2x/log(x) for x ≥ 55
    selberg_bound = 2 * x_test / np.log(x_test)
    print(f"  Selberg upper bound: {selberg_bound:.2f}")

    # Brun's sieve lower bound
    euler = np.prod([1 - 1/p for p in primerange(2, int(x_test**0.5) + 1)])
    print(f"  Euler product density: {euler:.4f}")
    print(f"  Expected coprime count: {x_test * euler:.2f}")

    # The gap between upper and lower bounds from sieve methods
    # is INHERENT (the "parity problem" of sieve theory).
    # Selberg's sieve cannot give an asymptotic better than 2x/log(x).
    # Bombieri's sieve gives π(x) ~ x/log(x) but with unspecified constant.

    print("\n  CONCLUSION: Sieve methods hit the PARITY BARRIER:")
    print("  they cannot distinguish a set from its complement.")
    print("  Exactly the same barrier that prevents polylog prime counting.")
    print("  The parity barrier IS the fundamental obstacle.")


def arithmetic_derivative():
    """
    WILD IDEA: The arithmetic derivative.

    Define n' for integers:
    - p' = 1 for prime p
    - (ab)' = a'b + ab' (Leibniz rule)
    - 0' = 0, 1' = 0

    Then p is prime iff p' = 1.

    n' has a closed form: n' = n * Σ_{p|n} v_p(n)/p
    where v_p(n) is the p-adic valuation.

    Computing n' requires knowing the factorization of n.

    But: n' = 1 iff n is prime.

    So: π(x) = |{n ≤ x : n' = 1}|

    Can we count the number of n ≤ x with n' = 1 without factoring each n?

    The set {n : n' = 1} is EXACTLY the set of primes.
    So this is a tautology.
    """
    print("\n" + "=" * 70)
    print("ARITHMETIC DERIVATIVE APPROACH")
    print("=" * 70)

    def arith_deriv(n):
        if n <= 1:
            return 0
        factors = factorint(n)
        return sum(n * e // p for p, e in factors.items())

    print("\nArithmetic derivatives:")
    for n in range(2, 31):
        nd = arith_deriv(n)
        print(f"  {n}' = {nd}", end="")
        if nd == 1:
            print(" ← PRIME")
        else:
            print()

    # Count n with n' = 1 up to x
    x = 100
    count = sum(1 for n in range(2, x+1) if arith_deriv(n) == 1)
    print(f"\n  |{{n ≤ {x} : n' = 1}}| = {count} = π({x}) = {primepi(x)}")

    print("\n  CONCLUSION: Arithmetic derivative encodes primality perfectly")
    print("  but computing it requires factorization — no shortcut.")


def number_field_sieve_analog():
    """
    IDEA: Counting primes in number fields.

    In Z[i] (Gaussian integers):
    - Primes p ≡ 1 (mod 4) split: p = π * π̄
    - Primes p ≡ 3 (mod 4) remain inert
    - p = 2 ramifies: 2 = -i(1+i)²

    The prime counting function for Z[i] involves BOTH rational primes
    and the arithmetic of the splitting. But it still relates to
    L-functions with zeros on the critical line.

    What about more exotic number fields?
    """
    print("\n" + "=" * 70)
    print("NUMBER FIELD APPROACH")
    print("=" * 70)

    # In Z[i]: prime ideals with norm ≤ x
    # Primes ≡ 1 (mod 4) contribute 2 ideals each
    # Primes ≡ 3 (mod 4) contribute 1 ideal each (with norm p²)
    # 2 contributes 1 ideal with norm 2

    x = 1000
    primes_mod4_1 = [p for p in primerange(2, x+1) if p % 4 == 1]
    primes_mod4_3 = [p for p in primerange(2, x+1) if p % 4 == 3]

    pi_x = primepi(x)
    pi_1 = len(primes_mod4_1)
    pi_3 = len(primes_mod4_3)

    print(f"\nPrime distribution mod 4 up to {x}:")
    print(f"  π(x) = {pi_x}")
    print(f"  p ≡ 1 (mod 4): {pi_1} ({pi_1/pi_x*100:.1f}%)")
    print(f"  p ≡ 3 (mod 4): {pi_3} ({pi_3/pi_x*100:.1f}%)")
    print(f"  p = 2: 1")

    # The "Chebyshev bias": typically more primes ≡ 3 (mod 4) than ≡ 1 (mod 4).
    print(f"  Bias (3 over 1): {pi_3 - pi_1}")

    # Key point: the splitting of primes in Z[i] is determined by
    # the Legendre symbol (-1/p), which is computable in O(log p) time.
    # But this tells us HOW a known prime splits, not WHERE primes are.

    # To count primes in Z[i] up to norm x, we still need π(x) for Q.
    # The "extra" information from the number field gives us the
    # mod-4 distribution, which is EASIER than counting all primes.

    # For general number field K of degree d:
    # π_K(x) = Li(x) + error involving zeros of ζ_K(s)
    # The error has the SAME structure (oscillatory, involves zeros).
    # More zeros (from the higher-degree L-functions) → harder, not easier.

    print("\n  CONCLUSION: Number fields DON'T make prime counting easier.")
    print("  The zeta zeros of K are a SUPERSET of those of Q.")
    print("  More algebraic structure = more zeros = more oscillations.")


if __name__ == "__main__":
    print("ALGEBRAIC/ARITHMETIC SHORTCUT EXPERIMENTS")
    print("Testing whether non-analytic approaches can count primes\n")

    wilson_based_counting()
    determinant_sieve()
    cyclotomic_approach()
    fast_factorial_mod()
    matrix_power_approach()
    quadratic_form_approach()
    arithmetic_derivative()
    number_field_sieve_analog()

    print("\n" + "=" * 70)
    print("GRAND SUMMARY")
    print("=" * 70)
    print("""
    Every algebraic approach tested reduces to one of:

    1. CIRCULAR: The formula requires knowing primes to compute primes.
       (Wilson, arithmetic derivative, cyclotomic)

    2. EQUIVALENT: The computation is equivalent to existing methods.
       (Determinant sieve ↔ Möbius inclusion-exclusion)

    3. HARDER: Adding algebraic structure makes the problem harder.
       (Number fields have MORE zeta zeros, not fewer)

    4. BARRIER: Fundamental obstacles prevent improvement.
       (Parity barrier in sieve theory, entropy of gap sequence)

    The common thread: PRIMALITY IS INHERENTLY LOCAL.
    There is no global algebraic identity that "solves" the primes.
    Every approach eventually requires checking individual numbers.

    The question becomes: Can we check FEWER numbers than O(x^{1/2})?
    """)
