"""
Session 14: Additive/multiplicative structure of primes for counting.

Key ideas to test:
1. Can we count primes by exploiting their distribution in residue classes?
   (Dirichlet: pi(x;q,a) ~ Li(x)/phi(q) for gcd(a,q)=1)
   If we could sum pi(x;q,a) for all admissible a, we get pi(x)-pi(q).

2. Can we use the Selberg sieve weights (quadratic optimization) to get
   a non-trivial POLYNOMIAL-TIME exact counting method?

3. Exploration of "smooth number" decomposition: every n <= x decomposes
   uniquely into smooth * rough parts. Can we count primes by counting
   specific rough numbers?

4. Can we express pi(x) as a sum over MULTIPLICATIVE characters?
   pi(x) = (1/phi(q)) sum_chi sum_{n<=x} chi(n) * [something]

5. NEW: "Batched primality" — if we could test O(polylog) batches of
   numbers simultaneously and combine results, could we count primes?
"""

import numpy as np
from sympy import primepi, isprime, primerange, factorint, totient
from math import isqrt, log, gcd
from collections import defaultdict


def residue_class_decomposition():
    """
    Decompose pi(x) into contributions from residue classes mod q.
    pi(x) = sum_{a coprime to q} pi(x; q, a) + pi(q)

    Question: is any individual pi(x; q, a) easier to compute than pi(x)?
    Or does their combination have cancelation that helps?
    """
    print("=" * 70)
    print("Residue class decomposition of pi(x)")
    print("=" * 70)

    for x in [100, 1000, 10000]:
        for q in [6, 30, 210]:
            # Count primes in each residue class
            counts = defaultdict(int)
            for p in primerange(2, x + 1):
                counts[p % q] += 1

            total = sum(counts.values())
            admissible = [a for a in range(q) if gcd(a, q) == 1]
            expected = total / len(admissible) if admissible else 0

            # Deviation from uniform
            deviations = {a: counts.get(a, 0) - expected for a in admissible}
            max_dev = max(abs(d) for d in deviations.values())
            rel_dev = max_dev / expected if expected > 0 else 0

            print(f"\nx={x}, q={q}: pi(x)={total}, {len(admissible)} admissible classes")
            print(f"  Expected per class: {expected:.1f}")
            print(f"  Max deviation: {max_dev:.1f} ({100*rel_dev:.1f}%)")
            print(f"  Chebyshev bias (1 vs q-1): {counts.get(1, 0)} vs {counts.get(q-1, 0)}")


def selberg_sieve_analysis():
    """
    The Selberg sieve gives upper bounds: pi(x) <= sum w_d * floor(x/d)
    where w_d are optimized quadratic weights.

    For the twin prime sieve: S = sum_{n<=x} (sum_{d|n(n+2)} lambda_d)^2
    where lambda_d are optimized to minimize S.

    Question: can Selberg weights give EXACT pi(x) (not just bounds)?
    The parity barrier says linear sieves can't distinguish even/odd
    number of factors. But QUADRATIC sieve weights are nonlinear.
    """
    print("\n" + "=" * 70)
    print("Selberg sieve weight analysis")
    print("=" * 70)

    for x in [100, 500, 1000]:
        sqrtx = isqrt(x)
        primes = list(primerange(2, sqrtx + 1))

        # Simple Selberg upper bound sieve
        # lambda_d = mu(d) * max(0, 1 - log(d)/log(z)) for squarefree d
        z = sqrtx
        log_z = log(z)

        # Compute sieve sum S = sum_{n<=x} (sum_{d|P(z), d|n} lambda_d)^2
        # which gives an upper bound for pi(x) - pi(z)

        # Direct computation: for each n, check if n is coprime to P(z)
        P_z = 1
        for p in primes:
            P_z *= p

        # Count unsieved numbers
        coprime_count = sum(1 for n in range(2, x + 1) if gcd(n, P_z) == 1 or n in primes)
        actual_pi = sum(1 for n in range(2, x + 1) if isprime(n))

        # The sieve remainder is: #{n <= x : gcd(n, P(z)) = 1} = pi(x) - pi(z) + 1
        unsieved = sum(1 for n in range(z + 1, x + 1) if all(n % p != 0 for p in primes))
        pi_z = len(primes)

        print(f"\nx={x}, z=sqrt(x)={sqrtx}")
        print(f"  pi(x) = {actual_pi}")
        print(f"  pi(z) = {pi_z}")
        print(f"  Unsieved in (z, x]: {unsieved}")
        print(f"  pi(x) - pi(z) + 1 = {actual_pi - pi_z + 1}")
        print(f"  Match: {unsieved == actual_pi - pi_z + 1}")

        # Key insight: the sieve EXACTLY counts primes (when z = sqrt(x))
        # because the unsieved numbers > sqrt(x) are exactly the primes
        # This is the Legendre sieve. But it costs 2^{pi(z)} terms.

        # Can we compute the unsieved count WITHOUT exponential inclusion-exclusion?
        # By Selberg: the upper bound is achieved by optimal weights
        # Selberg upper sieve gives: S <= 2*x/log(z) * (1 + o(1))
        selberg_bound = 2 * x / log_z if log_z > 0 else x
        print(f"  Selberg upper bound: {selberg_bound:.1f}")
        print(f"  Ratio actual/Selberg: {(actual_pi - pi_z)/selberg_bound:.3f}")


def smooth_rough_decomposition():
    """
    Every integer n has a unique decomposition n = s * r where s is the
    largest B-smooth part (all prime factors <= B) and r is the B-rough part
    (smallest prime factor > B).

    pi(x) counts the 1-rough numbers in [2, x] (i.e., primes = numbers with
    smallest factor = themselves, which is > 1).

    Can we count primes by first counting B-smooth numbers (easier) and
    then using the decomposition?
    """
    print("\n" + "=" * 70)
    print("Smooth/rough decomposition analysis")
    print("=" * 70)

    for x in [100, 500, 1000]:
        for B in [2, 3, 5, 10]:
            # Count B-rough primes (primes > B)
            rough_primes = sum(1 for p in primerange(B + 1, x + 1))

            # Count B-smooth numbers <= x
            # A number is B-smooth if all prime factors <= B
            smooth_count = 0
            for n in range(1, x + 1):
                factors = factorint(n)
                if all(p <= B for p in factors):
                    smooth_count += 1

            # Count numbers with smallest prime factor > B
            rough_count = 0
            for n in range(2, x + 1):
                smallest_factor = min(factorint(n).keys())
                if smallest_factor > B:
                    rough_count += 1

            print(f"\nx={x}, B={B}:")
            print(f"  B-smooth numbers: {smooth_count}")
            print(f"  B-rough numbers (smallest factor > B): {rough_count}")
            print(f"  Primes > B: {rough_primes}")
            print(f"  Rough but not prime: {rough_count - rough_primes}")

            # The "rough but not prime" count = numbers with all factors > B
            # but at least 2 factors. These are products of primes > B.
            # This is pi_2(x, B) + pi_3(x, B) + ... where pi_k counts k-almost-primes
            # with all factors > B.

            # For B = sqrt(x), rough_count = pi(x) - pi(B) + 1 (just primes > B and 1)
            # Actually rough_count also includes products p*q with p,q > B... NO
            # For B = sqrt(x), any rough number n <= x with smallest factor > sqrt(x)
            # must be prime (since n >= p^2 > x if n had two factors > sqrt(x)).
            if B >= isqrt(x):
                print(f"  B >= sqrt(x): rough = primes > B + {rough_count - rough_primes} correction")


def batched_primality_idea():
    """
    Can we test primality of MANY numbers simultaneously?

    Idea: Miller-Rabin with base 2 tests n by computing 2^{n-1} mod n.
    For a SET of numbers {n_1, ..., n_k}, can we batch-compute all
    2^{n_i - 1} mod n_i simultaneously faster than k individual tests?

    In TC^0, each test is independent and parallelizable. So testing k
    numbers takes the same depth as testing 1. The question is CIRCUIT SIZE:
    testing x numbers takes size O(x * poly(N)), which is 2^N * poly(N).

    Can we do better? E.g., is there shared structure in the computations
    2^{n-1} mod n for consecutive n?
    """
    print("\n" + "=" * 70)
    print("Batched primality / shared computation analysis")
    print("=" * 70)

    # For consecutive n, compute 2^{n-1} mod n
    # Note: 2^{n-1} mod n = 2^{n-1 mod lambda(n)} mod n where lambda is Carmichael
    # So the exponents are related to Carmichael function values

    x = 200
    results = []
    for n in range(2, x + 1):
        val = pow(2, n - 1, n)
        results.append((n, val, val == 1 or n == 2, isprime(n)))

    # How many composites pass MR base 2?
    pseudoprimes = [n for n, val, passes, is_p in results if passes and not is_p]
    print(f"MR base 2 pseudoprimes below {x}: {pseudoprimes}")

    # Analyze the pattern of 2^{n-1} mod n
    print(f"\n2^(n-1) mod n for n=2..30:")
    for n, val, passes, is_p in results[:29]:
        status = "PRIME" if is_p else ("PSP" if passes else f"={val}")
        print(f"  n={n:3d}: 2^{n-1} mod {n} = {val:4d}  [{status}]")

    # Key question: is there a way to COMBINE the MR tests for many n
    # without computing each independently?
    # E.g., can we compute sum_{n<=x} [2^{n-1} = 1 mod n] efficiently?

    # This sum counts primes + pseudoprimes. If we could compute it
    # in polylog time, and separately count pseudoprimes (which are rare),
    # we'd get pi(x).

    # Count: how many n <= x satisfy 2^{n-1} = 1 mod n?
    fermat_count = sum(1 for n, val, _, _ in results if val == 1)
    actual_primes = sum(1 for _, _, _, is_p in results if is_p)
    print(f"\n#{'{n<=x: 2^(n-1)=1 mod n}'} = {fermat_count}")
    print(f"pi({x}) = {actual_primes}")
    print(f"Pseudoprimes (Fermat): {fermat_count - actual_primes}")

    # The sum sum_{n<=x} [a^{n-1} = 1 mod n] is related to the
    # counting function for Fermat pseudoprimes. Are these easier to count?
    # Erdos (1950): the number of pseudoprimes to base 2 below x is
    # x^{1-c/ln ln x} for some c > 0. So they're rare but infinitely many.

    # For computing pi(x), this approach would need:
    # 1. Efficiently compute #{n<=x : 2^{n-1}=1 mod n} (the Fermat count)
    # 2. Efficiently compute #{pseudoprimes <= x}
    # Then pi(x) = Fermat count - pseudoprime count

    # Problem: step 1 is as hard as pi(x) because it requires testing each n.
    # No known way to batch-evaluate 2^{n-1} mod n for all n simultaneously
    # in less than O(x) operations.
    print("\nConclusion: batched Fermat/MR is O(x * poly(N)) — no savings over pi(x)")


def character_sum_approach():
    """
    Can we express pi(x) via Dirichlet characters?

    pi(x) = sum_{n<=x} 1_prime(n)

    For n > q, n coprime to q:
    1_{n in residue class a mod q} = (1/phi(q)) sum_chi chi(a_bar) chi(n)

    So: pi(x) = pi(q) + sum_{a coprime q} (1/phi(q)) sum_chi chi(a_bar) sum_{n<=x} chi(n) * 1_prime(n)
    = pi(q) + (1/phi(q)) sum_chi [sum_{p<=x} chi(p)] * [sum_a chi(a_bar)]

    But sum_a chi(a_bar) = phi(q) if chi = chi_0 (principal), 0 otherwise.
    So this just gives pi(x) = pi(q) + pi(x) - pi(q). Tautological.

    What about using the EXPLICIT FORMULA for each L(s, chi)?
    sum_{p<=x} chi(p) involves zeros of L(s, chi).
    """
    print("\n" + "=" * 70)
    print("Character sum approach")
    print("=" * 70)

    # Test: for various q, how well can we predict pi(x) from the
    # distribution across residue classes?
    x = 1000
    for q in [6, 30, 210, 2310]:
        residues = defaultdict(int)
        for p in primerange(q + 1, x + 1):
            residues[p % q] += 1

        admissible = sorted([a for a in range(q) if gcd(a, q) == 1])
        n_admissible = len(admissible)

        # If primes were equidistributed, each class would have pi(x)/phi(q) primes
        expected = (int(primepi(x)) - int(primepi(q))) / n_admissible

        # Compute variance of deviations
        deviations = [residues.get(a, 0) - expected for a in admissible]
        variance = np.var(deviations)
        max_dev = max(abs(d) for d in deviations)

        # Chebyshev bias: count pi(x;q,NR) vs pi(x;q,R) where NR/R are
        # non-residues vs residues modulo q
        print(f"\nq={q}, phi(q)={n_admissible}")
        print(f"  Expected per class: {expected:.2f}")
        print(f"  Variance of deviations: {variance:.2f}")
        print(f"  Max deviation: {max_dev:.1f} ({100*max_dev/expected:.1f}%)")
        print(f"  Info gain over uniform: {variance/expected:.4f} (should → 0)")

    print("\nConclusion: character sums give asymptotic equidistribution, not exact counts")
    print("Each L(s,chi) has its own zeros, making the problem HARDER not easier")


if __name__ == '__main__':
    residue_class_decomposition()
    selberg_sieve_analysis()
    smooth_rough_decomposition()
    batched_primality_idea()
    character_sum_approach()
